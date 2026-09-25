import argparse
import pathlib
import random
import re
import string
from unittest.mock import patch

import pytest

from AletheiaSeq import (
    aletheiaseq_commandline,  # import argument_parser, _is_valid_path, _check_folder, _is_valid_blast
    aletheiaseq_parse,
    aletheiaseq_prepare,
)

"""
FIXTURES
"""


@pytest.fixture
def temp_fasta_file(tmp_path):
    file_path = tmp_path / "test_seqs.fasta"
    file_path.write_text(
        ">seq1\nACTGCACTGACTTCCAGTCGACCGAGCTTT\n>seq2\nACGAAACTATCCGAGCTTAGCTAGCTTTAGCGGGATCGCGCGAGACT"
    )
    return str(file_path)


@pytest.fixture
def blast_db(tmp_path):
    for ext in [".nhr", ".nin", ".nsq", ".njs"]:
        (tmp_path / f"db{ext}").touch()

    return str(tmp_path / "db")


@pytest.fixture
def valid_yaml(tmp_path):
    file_path = tmp_path / "test.yaml"
    file_path.touch()

    return str(file_path)


@pytest.fixture
def base_args(valid_yaml, tmp_path):
    rdm_str = "".join(random.choice(string.ascii_letters) for _ in range(12))
    return [
        "--yaml",
        valid_yaml,
        "--outfolder",
        str(tmp_path / rdm_str),
    ]


@pytest.fixture
def blast_out(tmp_path):
    blast_path = tmp_path / "blast_out.tsv"
    blast_path.write_text("qseqid,sseqid,qlen,slen,qstart,qend,sstart,send,evalue,bitscore,length,pident\n")

    return str(blast_path)


"""
tests for argparse type checks
"""


def test_is_valid_path_pass(temp_fasta_file):
    assert aletheiaseq_commandline._is_valid_path(str(temp_fasta_file)) == str(temp_fasta_file)


def test_is_valid_path_fail():
    with pytest.raises(argparse.ArgumentTypeError, match="is not a valid path"):
        aletheiaseq_commandline._is_valid_path("/not/a/valid/path.fasta")


def test_check_folder_create(tmp_path):
    new_dir = tmp_path / "outfolder"
    assert aletheiaseq_commandline._check_folder(str(new_dir)) == str(new_dir)
    assert new_dir.is_dir()


def test_check_folder_existing_dir(tmp_path):
    with pytest.warns(UserWarning, match="Output directory already exists"):
        assert aletheiaseq_commandline._check_folder(str(tmp_path)) == str(tmp_path)


def test_check_folder_fail(temp_fasta_file):
    with pytest.raises(argparse.ArgumentTypeError, match="is a path not a directory"):
        aletheiaseq_commandline._check_folder(str(temp_fasta_file))


def test_check_folder_permission_error(tmp_path):
    with (
        patch.object(pathlib.Path, "mkdir", side_effect=PermissionError("Permission denied")),
        pytest.raises(argparse.ArgumentTypeError, match="invalid or a system I/O error occurred"),
    ):
        aletheiaseq_commandline._check_folder(str(tmp_path / "outfolder"))


def test_check_folder_missing_parent(tmp_path):
    wrong_path = tmp_path / "parent_does_not_exist" / "outfolder"
    with pytest.raises(argparse.ArgumentTypeError, match="missing required parent folder"):
        aletheiaseq_commandline._check_folder(str(wrong_path))


def test_valid_blast(blast_db):
    assert aletheiaseq_commandline._is_valid_blast(blast_db) == blast_db


@pytest.mark.parametrize("missing_ext", [".nhr", ".nin", ".nsq", ".njs"])
def test_missing_blast_file(tmp_path, missing_ext):
    for ext in [".nhr", ".nin", ".nsq", ".njs"]:
        if ext != missing_ext:
            (tmp_path / f"db{ext}").touch()

    with pytest.raises(
        argparse.ArgumentTypeError, match=rf"BlastDB:{re.escape(missing_ext)} is not present in directory given"
    ):
        aletheiaseq_commandline._is_valid_blast(str(tmp_path / "db"))


"""
tests for argparse
"""


def test_argument_parser_constructs():
    parser = aletheiaseq_commandline.argument_parser()
    assert parser.prog == "aletheiaseq"


@pytest.mark.parametrize("missing_arg", ["--yaml", "--outfolder"])
def test_requires_top_level_args(missing_arg, valid_yaml, tmp_path, blast_out, capsys):
    arg_pairs = [["--yaml", valid_yaml], ["--outfolder", str(tmp_path / missing_arg)]]

    args = []
    for pair in arg_pairs:
        if pair[0] != missing_arg:
            args += pair

    args += ["parse", "--blast_out", blast_out, "--sample_id", "ID-123"]

    with pytest.raises(SystemExit) as exc:
        aletheiaseq_commandline.argument_parser().parse_args(args)

    captured = capsys.readouterr()
    assert "required" in captured.err
    assert missing_arg in captured.err
    assert exc.value.code == 2


def test_subcommand_required(base_args, capsys):
    parser = aletheiaseq_commandline.argument_parser()

    with pytest.raises(SystemExit) as exc:
        parser.parse_args(base_args)

    captured = capsys.readouterr()
    assert "required" in captured.err
    assert "subcommand" in captured.err
    assert exc.value.code == 2

    with pytest.raises(SystemExit) as exc, pytest.warns(UserWarning, match="Output directory already exists"):
        parser.parse_args(base_args + ["prepare"])

    captured = capsys.readouterr()
    assert "required" in captured.err
    assert "--fasta" in captured.err
    assert exc.value.code == 2

    with pytest.raises(SystemExit) as exc, pytest.warns(UserWarning, match="Output directory already exists"):
        parser.parse_args(base_args + ["parse"])

    captured = capsys.readouterr()
    assert "required" in captured.err
    assert "--sample_id" in captured.err
    assert exc.value.code == 2


def test_prepare_subcommand_parses(base_args, temp_fasta_file, blast_db):
    parser = aletheiaseq_commandline.argument_parser()
    args = parser.parse_args(
        base_args
        + [
            "prepare",
            "--fasta",
            temp_fasta_file,
            "--blast_db",
            blast_db,
        ]
    )

    assert args.subcommand == "prepare"
    assert args.yaml == base_args[1]
    assert args.blast_db == blast_db
    assert args.out_folder == base_args[3]
    assert args.ref_fasta == temp_fasta_file
    assert args.skope_idx is None


def test_prepare_subcommand_optional(base_args, temp_fasta_file, tmp_path, blast_db):
    parser = aletheiaseq_commandline.argument_parser()
    (tmp_path / "skope.idx").touch()
    args = parser.parse_args(
        base_args
        + [
            "prepare",
            "--fasta",
            temp_fasta_file,
            "--blast_db",
            blast_db,
        ]
        + ["--skope_idx", str(tmp_path / "skope.idx")]
    )

    assert args.skope_idx == str(tmp_path / "skope.idx")


def test_prepare_subcommand_function(base_args, temp_fasta_file, blast_db):
    parser = aletheiaseq_commandline.argument_parser()
    args = parser.parse_args(
        base_args
        + [
            "prepare",
            "--fasta",
            temp_fasta_file,
            "--blast_db",
            blast_db,
        ]
    )

    assert args.func is aletheiaseq_prepare.prepare


def test_parse_subcommand_parses(base_args, blast_out):
    parser = aletheiaseq_commandline.argument_parser()
    args = parser.parse_args(
        base_args
        + [
            "parse",
            "--sample_id",
            "ID-123",
            "--blast_out",
            blast_out,
        ]
    )

    assert args.subcommand == "parse"
    assert args.sample == "ID-123"
    assert args.blast_out == blast_out
    assert args.skope_out is None
    assert args.onyx_server is None


def test_parse_subcommand_optional(base_args, tmp_path, blast_out):
    parser = aletheiaseq_commandline.argument_parser()
    (tmp_path / "skope_classify.out").touch()
    args = parser.parse_args(
        base_args
        + [
            "parse",
            "--sample_id",
            "ID-123",
            "--blast_out",
            blast_out,
        ]
        + [
            "--skope_out",
            str(tmp_path / "skope_classify.out"),
        ]
        + ["--onyx_server", "mscape"]
    )

    assert args.skope_out == str(tmp_path / "skope_classify.out")
    assert args.onyx_server == "mscape"


def test_parse_subcommand_function(base_args, blast_out):
    parser = aletheiaseq_commandline.argument_parser()
    args = parser.parse_args(
        base_args
        + [
            "parse",
            "--sample_id",
            "ID-123",
            "--blast_out",
            blast_out,
        ]
    )

    assert args.func is aletheiaseq_parse.parse


def test_publish_flag(base_args, blast_out):
    parser = aletheiaseq_commandline.argument_parser()
    args = parser.parse_args(base_args + ["parse", "--sample_id", "ID-123", "--blast_out", blast_out, "--publish"])

    assert args.publish
    assert args.dryrun


def test_dryrun_flag(base_args, blast_out):
    parser = aletheiaseq_commandline.argument_parser()
    args = parser.parse_args(
        base_args + ["parse", "--sample_id", "ID-123", "--blast_out", blast_out, "--publish", "--no_dryrun"]
    )

    assert args.publish
    assert not args.dryrun
