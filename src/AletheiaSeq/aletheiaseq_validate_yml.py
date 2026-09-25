import ast
import logging
import re
import warnings
from datetime import date
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, Field, ValidationError, field_validator, model_validator


class SchemaDetails(BaseModel):
    name: str
    schema_type: Literal["speciation", "detection"]
    version: str

    @field_validator("name")
    @classmethod
    def validate_name(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("Schema name cannot be empty")

        return value

    @field_validator("version")
    @classmethod
    def check_semantic_versioning(cls, value: str) -> str:
        if not re.match(r"^\d+\.\d+\.\d+$", value):
            message = f"WARN: version ({value}) is not in the semantic versioning format <major>.<minor>.<patch>. This will not prevent AletheiaSeq from running but is good practice."
            warnings.warn(
                message,
                UserWarning,
                stacklevel=2,
            )

        return value


class BlastDetails(BaseModel):
    seqs: list[str] = Field(
        min_length=1,
        description="Sequences contained within the blast database",
    )
    db_name: str
    db_type: Literal["nucleotide", "protein"]
    last_update: date
    outfmt: list[str]
    filters: list[str]
    min_hits: int = Field(
        default=3,
        ge=1,
        description="Number of times loci must be seen in the input data to be considered a detection",
    )

    @field_validator("seqs")
    @classmethod
    def validate_unique_seqs(cls, value: list[str]) -> list[str]:
        if len(value) != len(set(value)):
            raise ValueError("Duplicate loci found in BlastDetails[seqs]")

        return value

    @field_validator("last_update")
    @classmethod
    def validate_last_update(cls, value: date) -> date:
        if value > date.today():
            raise ValueError("BlastDetails[last_update] cannot be in the future")

        return value

    @model_validator(mode="after")
    def validate_allowed_outfmt_fields(self):
        allowed_columns = {
            "qseqid",
            "sseqid",
            "qlen",
            "slen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "length",
            "pident",
        }

        incorrect = set(self.outfmt) - allowed_columns

        if incorrect:
            raise ValueError(f"BlastDetails[outfmt] contains invalid column(s): {sorted(incorrect)}")

        return self

    @model_validator(mode="after")
    def validate_required_outfmt_fields(self):
        required_columns = {"qseqid", "sseqid", "slen", "evalue", "bitscore", "length", "pident"}

        missing = required_columns - set(self.outfmt)

        if missing:
            raise ValueError(f"BlastDetails[outfmt] missing required column(s): {sorted(missing)}")

        return self

    @field_validator("outfmt")
    @classmethod
    def validate_unique_outfmt(cls, value: list[str]) -> list:
        if len(value) != len(set(value)):
            raise ValueError("Duplicate columns in BlastDetails[outfmt]")

        return value

    def process_equation_side(self, node: ast.AST) -> bool:
        # check if valid column name
        if isinstance(node, ast.Name):
            if node.id not in self.outfmt:
                raise ValueError(f"Column used to filter is not contained in blast outfmt list: {node.id}")
            else:
                return True
        # check if numeric value
        elif isinstance(node, ast.Constant):
            if not isinstance(node.value, int | float) or isinstance(node.value, bool):
                raise ValueError(f"Only numeric constants are permitted in blast filter: {ast.unparse(node)}")
            else:
                return False
        elif isinstance(node, ast.UnaryOp):
            if not isinstance(node.op, ast.UAdd | ast.USub):
                raise ValueError(f"Unsupported unary operator in blast filter: {ast.unparse(node)}")
            else:
                # process the non +/- bit of the equation to make sure it is a valid number or column name
                return self.process_equation_side(node.operand)
        # check if arithmetic expression e.g., (col_name * 0.8)
        elif isinstance(node, ast.BinOp):
            allowed_operators = (
                ast.Add,
                ast.Sub,
                ast.Mult,
                ast.Div,
            )
            if not isinstance(node.op, allowed_operators):
                raise ValueError(f"Unsupported arithmetic operator in blast filter: {ast.unparse(node)}")
            left_has_column = self.process_equation_side(node.left)
            right_has_column = self.process_equation_side(node.right)

            return left_has_column or right_has_column
        else:
            raise ValueError(f"Expression not supported for blast results filter: {ast.unparse(node)}")

    @model_validator(mode="after")
    def validate_filter_expressions(self):
        filter_list = self.filters
        parsed_filters = []
        for filt in filter_list:
            try:
                parsed_equation = ast.parse(filt, mode="eval")
            except SyntaxError as exc:
                raise ValueError(f"Invalid query syntax in baslt filter: {filt}") from exc
            for node in ast.walk(parsed_equation):
                if isinstance(node, ast.Call | ast.Attribute):
                    raise ValueError(f"Function calls and Attribute access are not supported: {ast.unparse(node)}")

            if isinstance(parsed_equation.body, ast.BinOp) and isinstance(
                parsed_equation.body.op, ast.BitOr | ast.BitAnd
            ):
                # this means there is | or & in the filter so both sides of the filter need to be processed
                parsed_filters += [parsed_equation.body.left, parsed_equation.body.right]
            else:
                parsed_filters.append(parsed_equation.body)

        for equation in parsed_filters:
            if not isinstance(equation, ast.Compare):
                raise ValueError(f"Expected a comparison to filter blast results: {ast.unparse(equation)}")
            elif len(equation.ops) > 1 or len(equation.comparators) > 1:
                raise ValueError(f"Chained comparisons are not supported: {ast.unparse(equation)}")
            else:
                allowed_operators = (
                    ast.Eq,
                    ast.NotEq,
                    ast.Lt,
                    ast.LtE,
                    ast.Gt,
                    ast.GtE,
                )

                operator = equation.ops[0]
                if not isinstance(operator, allowed_operators):
                    raise ValueError(
                        f"Unsupported operator in blast filter: {type(operator).__name__}. "
                        "Allowed operators are ==, !=, <, <=, >, >="
                    )

                left = equation.left
                right = equation.comparators[0]

                at_least_one_column_left = self.process_equation_side(left)
                at_least_one_column_right = self.process_equation_side(right)

                if not (at_least_one_column_left or at_least_one_column_right):
                    raise ValueError(f"No valid column name is present in blast filter: {filt}")

        return self


class SchemaGroup(BaseModel):
    name: str
    defining_loci: list[str]
    required_loci: list[str]
    minimum_loci_required: int = Field(
        default=1,
        ge=1,
        description="Number of loci required to be deceted above QC thresholds to return positive result for group",
    )

    @field_validator("name")
    @classmethod
    def validate_schemagroup_name(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("SchemaGroup name cannot be empty")

        return value

    @field_validator("defining_loci")
    @classmethod
    def validate_seqs(cls, value: list[str]) -> list[str]:
        if not value or len(value) == 0:
            raise ValueError("At least one locus is required in SchemaGroup[defining_loci]")

        return value

    @field_validator("defining_loci")
    @classmethod
    def validate_unique_defining_seqs(cls, value: list[str]) -> list[str]:
        if len(value) != len(set(value)):
            raise ValueError("Duplicate loci found in SchemaGroups[defining_loci]")

        return value

    @field_validator("required_loci")
    @classmethod
    def validate_unique_required_seqs(cls, value: list[str]) -> list[str]:
        if len(value) != len(set(value)):
            raise ValueError("Duplicate loci found in SchemaGroups[required_loci]")

        return value

    @model_validator(mode="after")
    def validate_required_subset(self):
        missing = set(self.required_loci) - set(self.defining_loci)

        if missing:
            raise ValueError(f"required_loci must be a subset of defining_loci. Invalid: {sorted(missing)}")

        return self

    @model_validator(mode="after")
    def validate_minimum(self):
        if self.minimum_loci_required > len(self.defining_loci):
            raise ValueError("minimum_loci_required exceeds number of defining_loci")

        if len(self.required_loci) > self.minimum_loci_required:
            raise ValueError("minimum_loci_required cannot be less than the number of required_loci")

        return self


class SchemaConfig(BaseModel):
    SchemaDetails: SchemaDetails
    BlastDetails: BlastDetails
    SchemaGroups: list[SchemaGroup] = Field(min_length=1)

    @model_validator(mode="after")
    def validate_group_loci_exist(self):
        blast_loci = set(self.BlastDetails.seqs)

        for group in self.SchemaGroups:
            missing = set(group.defining_loci) - blast_loci

            if missing:
                raise ValueError(
                    f"Group '{group.name}' references loci not listed in BlastDetails.seqs: {sorted(missing)}"
                )

        return self

    @model_validator(mode="after")
    def validate_group_names_unique(self):
        names = [g.name.casefold() for g in self.SchemaGroups]

        if len(names) != len(set(names)):
            raise ValueError("SchemaGroup names must be unique")

        return self

    @model_validator(mode="after")
    def validate_unused_loci(self):
        used = set()

        for group in self.SchemaGroups:
            used.update(group.defining_loci)

        unused = set(self.BlastDetails.seqs) - used

        if unused:
            raise ValueError(f"Unused loci in BlastDetails.seqs: {sorted(unused)}")

        return self

    def summarise_methods(self):
        blast_filters = {
            "filters": self.BlastDetails.filters,
            "minimum_blast_hits": self.BlastDetails.min_hits,
        }
        groups = {}
        for g in self.SchemaGroups:
            groups[g.name] = g.model_dump(exclude=["name"])

        return {"blast_filters": blast_filters, "schema_group_definitions": groups}


def _load_yml(yaml_file: str):
    with Path(yaml_file).open("r") as f:
        y = yaml.safe_load(f)

    return y


def validate_yaml(yaml_path: str) -> SchemaConfig:
    yml = _load_yml(yaml_path)
    logging.info("Yaml file loaded, validating contents.")

    try:
        config = SchemaConfig.model_validate(yml)
    except ValidationError as e:
        logging.info(f"Schema: {config.SchemaDetails.name}, v{config.SchemaDetails.version} validation failed.")
        logging.info(e.errors())

    logging.info(f"Schema: {config.SchemaDetails.name}, v{config.SchemaDetails.version} validated.")
    return config
