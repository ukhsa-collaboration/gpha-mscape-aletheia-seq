import json
from datetime import date, timedelta

import pytest

from AletheiaSeq import aletheiaseq_validate_yml


@pytest.fixture
def valid_schema_details():
    return {
        "name": "schema1",
        "schema_type": "speciation",
        "version": "0.0.1",
    }


@pytest.fixture
def valid_blast_details():
    return {
        "seqs": ["locus1", "locus2"],
        "db_name": "test_db",
        "db_type": "nucleotide",
        "last_update": date.today(),
        "outfmt": [
            "qseqid",
            "sseqid",
            "qlen",
            "slen",
            "evalue",
            "bitscore",
            "length",
            "pident",
        ],
        "filters": ["pident >= 80"],
        "min_hits": 3,
    }


def build_blast_details(config, filt):
    config = config.copy()
    config["filters"] = [filt]
    return aletheiaseq_validate_yml.BlastDetails.model_validate(config)


def build_blast_seqs(config, seqs):
    config = config.copy()
    config["seqs"] = seqs
    return aletheiaseq_validate_yml.BlastDetails.model_validate(config)


@pytest.fixture
def valid_schema_group():
    return {
        "name": "schema group 1",
        "defining_loci": ["locus1", "locus2"],
        "required_loci": ["locus1"],
        "minimum_loci_required": 1,
    }


@pytest.fixture
def valid_schema_config(valid_schema_details, valid_blast_details, valid_schema_group):
    return {
        "SchemaDetails": valid_schema_details,
        "BlastDetails": valid_blast_details,
        "SchemaGroups": [valid_schema_group],
    }


@pytest.fixture
def temp_yml_file(tmp_path, valid_schema_config):
    file_path = tmp_path / "test_schema.yml"
    file_path.write_text(
        json.dumps(
            valid_schema_config, default=str
        )  # this is to convert the datetime in blast details and is a bit of a hacky way of doing it, soz
    )
    return str(file_path)


"""
TEST SchemaDetails
"""


def test_schema_details_validate(valid_schema_details):
    validated_model = aletheiaseq_validate_yml.SchemaDetails.model_validate(valid_schema_details)
    assert isinstance(validated_model, aletheiaseq_validate_yml.SchemaDetails)


def test_incorrect_schema_type(valid_schema_details):
    config = valid_schema_details.copy()
    config["schema_type"] = "random string"

    with pytest.raises(
        ValueError,
    ):
        aletheiaseq_validate_yml.SchemaDetails.model_validate(config)


def test_empty_schema_name(valid_schema_details):
    config = valid_schema_details.copy()
    config["name"] = ""

    with pytest.raises(
        ValueError,
    ):
        aletheiaseq_validate_yml.SchemaDetails.model_validate(config)


def test_warn_sematic_versioning(valid_schema_details):
    config = valid_schema_details.copy()
    config["version"] = "version 1"

    with pytest.warns(UserWarning, match="semantic versioning format"):
        aletheiaseq_validate_yml.SchemaDetails.model_validate(config)


"""
TEST BlastDetails
"""


@pytest.mark.parametrize(
    "seqs",
    [[], ["locus1", "locus1"]],
)
def test_seqs_check_fail(valid_blast_details, seqs):
    with pytest.raises(
        ValueError,
    ):
        build_blast_seqs(valid_blast_details, seqs)


def test_invalid_db_type(valid_blast_details):
    config = valid_blast_details.copy()
    config["db_type"] = "invalid db type"
    with pytest.raises(
        ValueError,
    ):
        aletheiaseq_validate_yml.BlastDetails.model_validate(config)


def test_default_min_hits(valid_blast_details):
    config = valid_blast_details.copy()
    del config["min_hits"]

    validated_model = aletheiaseq_validate_yml.BlastDetails.model_validate(config)
    assert validated_model.min_hits == 3


def test_minimum_min_hits(valid_blast_details):
    config = valid_blast_details.copy()
    config["min_hits"] = 0

    with pytest.raises(
        ValueError,
    ):
        aletheiaseq_validate_yml.BlastDetails.model_validate(config)


def test_last_update_future_date(valid_blast_details):
    config = valid_blast_details.copy()
    config["last_update"] += timedelta(days=7)

    with pytest.raises(ValueError, match="cannot be in the future"):
        aletheiaseq_validate_yml.BlastDetails.model_validate(config)


@pytest.mark.parametrize(
    "outfmt",
    [
        [  # missing required column
            "qseqid",
            "sseqid",
            "qlen",
            "slen",
            "evalue",
            "bitscore",
            "length",
        ],
        [  # includes unknown column name
            "qseqid",
            "sseqid",
            "qlen",
            "slen",
            "evalue",
            "bitscore",
            "length",
            "pident",
            "unknown column",
        ],
        [  # includes duplicated column name
            "qseqid",
            "sseqid",
            "qlen",
            "slen",
            "evalue",
            "bitscore",
            "length",
            "pident",
            "qstart",
            "qstart",
        ],
    ],
)
def test_outfmt_column_requirements(valid_blast_details, outfmt):
    config = valid_blast_details.copy()
    config["outfmt"] = outfmt

    with pytest.raises(
        ValueError,
    ):
        aletheiaseq_validate_yml.BlastDetails.model_validate(config)


"""
TEST BlastDetails filtering
"""

"""
Should pass
"""


@pytest.mark.parametrize(
    "filter_expression",
    [
        "pident >= 80",
        "length <= slen",
        "length >= (slen * 0.8)",
        "length >= ((slen + qlen) / 2)",
        "evalue <= 1e-10",
        "bitscore >= -10",
        "(length <= slen) | (length >= (0.8 * slen))",
    ],
)
def test_valid_filter_expressions(valid_blast_details, filter_expression):
    result = build_blast_details(
        valid_blast_details,
        filter_expression,
    )

    assert result.filters == [filter_expression]


"""
Should fail
"""


def test_filter_unknown_column(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Column used to filter",
    ):
        build_blast_details(
            valid_blast_details,
            "unknown_column >= 80",
        )


def test_filter_no_column_present(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="No valid column name",
    ):
        build_blast_details(
            valid_blast_details,
            "80 >= 70",
        )


def test_filter_chained_comparison(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Chained comparisons",
    ):
        build_blast_details(
            valid_blast_details,
            "pident > qlen > 80",
        )


def test_filter_unsupported_operator(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Unsupported operator",
    ):
        build_blast_details(
            valid_blast_details,
            "pident in [80, 90]",
        )


def test_filter_function_call_not_allowed(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Function calls",
    ):
        build_blast_details(
            valid_blast_details,
            "len(pident) > 10",
        )


def test_filter_attribute_access_not_allowed(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Attribute access",
    ):
        build_blast_details(
            valid_blast_details,
            "obj.attribute > 10",
        )


def test_filter_invalid_syntax(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Invalid query syntax",
    ):
        build_blast_details(
            valid_blast_details,
            "pident >=",
        )


def test_filter_invalid_arithmetic_operator(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Unsupported arithmetic operator",
    ):
        build_blast_details(
            valid_blast_details,
            "length >= (slen ** 2)",
        )


def test_filter_string_constant_not_allowed(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Only numeric constants",
    ):
        build_blast_details(
            valid_blast_details,
            "pident >= 'high'",
        )


@pytest.mark.parametrize(
    "filter_expression",
    [
        "pident + 80",
        "slen * 0.8",
        "42",
        "'test'",
        "[1, 2, 3]",
    ],
)
def test_filter_must_be_comparison(valid_blast_details, filter_expression):
    with pytest.raises(
        ValueError,
        match="Expected a comparison to filter blast results",
    ):
        build_blast_details(
            valid_blast_details,
            filter_expression,
        )


def test_filter_boolean_expression_not_supported(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Expected a comparison to filter blast results",
    ):
        build_blast_details(
            valid_blast_details,
            "pident >= 80 and length >= 100",
        )


def test_filter_rejects_unsupported_expression_type(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Expression not supported for blast results filter",
    ):
        build_blast_details(
            valid_blast_details,
            "pident >= [80]",
        )


def test_filter_rejects_bitwise_not_operator(valid_blast_details):
    with pytest.raises(
        ValueError,
        match="Unsupported unary operator",
    ):
        build_blast_details(
            valid_blast_details,
            "pident >= ~1",
        )


"""
TEST SchemaGroup
"""


def test_name_not_empty(valid_schema_group):
    config = valid_schema_group.copy()
    config["name"] = ""

    with pytest.raises(ValueError, match="cannot be empty"):
        aletheiaseq_validate_yml.SchemaGroup.model_validate(config)


def test_min_loci_default(valid_schema_group):
    config = valid_schema_group.copy()
    del config["minimum_loci_required"]

    validated_model = aletheiaseq_validate_yml.SchemaGroup.model_validate(config)
    assert validated_model.minimum_loci_required == 1


def test_min_loci_minimum(valid_schema_group):
    config = valid_schema_group.copy()
    config["minimum_loci_required"] = 0

    with pytest.raises(
        ValueError,
    ):
        aletheiaseq_validate_yml.SchemaGroup.model_validate(config)


def test_required_not_in_defining(valid_schema_group):
    config = valid_schema_group.copy()
    config["defining_loci"] = [
        "locus2",
    ]

    with pytest.raises(ValueError, match="must be a subset"):
        aletheiaseq_validate_yml.SchemaGroup.model_validate(config)


def test_check_empty_loci(valid_schema_group):
    config = valid_schema_group.copy()
    config["defining_loci"] = []
    config["required_loci"] = []

    with pytest.raises(ValueError, match="one locus is required"):
        aletheiaseq_validate_yml.SchemaGroup.model_validate(config)


def test_check_duplicate_defining_loci(valid_schema_group):
    config = valid_schema_group.copy()
    config["defining_loci"] = [
        "locus1",
        "locus1",
    ]
    config["required_loci"] = []

    with pytest.raises(ValueError, match="Duplicate loci found"):
        aletheiaseq_validate_yml.SchemaGroup.model_validate(config)


def test_check_duplicate_required_loci(valid_schema_group):
    config = valid_schema_group.copy()
    config["required_loci"] = [
        "locus1",
        "locus1",
    ]

    with pytest.raises(ValueError, match="Duplicate loci found"):
        aletheiaseq_validate_yml.SchemaGroup.model_validate(config)


def test_min_loci_larger_than_defining(valid_schema_group):
    config = valid_schema_group.copy()
    config["minimum_loci_required"] = 3

    with pytest.raises(ValueError, match="minimum_loci_required exceeds"):
        aletheiaseq_validate_yml.SchemaGroup.model_validate(config)


def test_min_loci_less_than_required(valid_schema_group):
    config = valid_schema_group.copy()
    config["required_loci"] = [
        "locus1",
        "locus2",
    ]

    with pytest.raises(ValueError, match="less than the number of required_loci"):
        aletheiaseq_validate_yml.SchemaGroup.model_validate(config)


"""
TEST SchemaConfig
"""


def test_valid_schema_config(valid_schema_config):
    valid_config_model = aletheiaseq_validate_yml.SchemaConfig.model_validate(valid_schema_config)

    assert valid_config_model.SchemaDetails == aletheiaseq_validate_yml.SchemaDetails.model_validate(
        valid_schema_config["SchemaDetails"]
    )
    assert valid_config_model.BlastDetails == aletheiaseq_validate_yml.BlastDetails.model_validate(
        valid_schema_config["BlastDetails"]
    )
    assert valid_config_model.SchemaGroups[0] == aletheiaseq_validate_yml.SchemaGroup.model_validate(
        valid_schema_config["SchemaGroups"][0]
    )


def test_no_Schema_groups(valid_schema_config):
    config = valid_schema_config.copy()
    config["SchemaGroups"] = []

    with pytest.raises(
        ValueError,
    ):
        aletheiaseq_validate_yml.SchemaConfig.model_validate(config)


def test_schema_group_loci_not_in_blast_details(valid_schema_config):
    config = valid_schema_config.copy()
    config["BlastDetails"]["seqs"] = ["locus1"]

    with pytest.raises(ValueError, match="references loci not listed"):
        aletheiaseq_validate_yml.SchemaConfig.model_validate(config)


def test_non_unique_schema_group(valid_schema_config, valid_schema_group):
    config = valid_schema_config.copy()
    config["SchemaGroups"].append(valid_schema_group)

    with pytest.raises(ValueError, match="names must be unique"):
        aletheiaseq_validate_yml.SchemaConfig.model_validate(config)


def test_unused_loci_in_blast_seqs(valid_schema_config):
    config = valid_schema_config.copy()
    config["BlastDetails"]["seqs"].append("locus3")

    with pytest.raises(ValueError, match="Unused loci in BlastDetails.seqs"):
        aletheiaseq_validate_yml.SchemaConfig.model_validate(config)


"""
TEST read file
"""


def test_read_from_file_and_validate(temp_yml_file, valid_schema_config):
    validated_model = aletheiaseq_validate_yml.validate_yaml(temp_yml_file)

    assert validated_model == aletheiaseq_validate_yml.SchemaConfig.model_validate(valid_schema_config)
