"""Packaged CLI entrypoint for DrugReflector."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .drug_reflector import DrugReflector
from .utils import load_h5ad_file


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="DrugReflector inference tools",
    )
    parser.add_argument(
        "input",
        nargs="?",
        help="Input H5AD file for batch prediction, or omit when using `serve`.",
    )
    parser.add_argument(
        "--model1",
        type=str,
        default="checkpoints/model_fold_0.pt",
        help="Path to first model checkpoint.",
    )
    parser.add_argument(
        "--model2",
        type=str,
        default="checkpoints/model_fold_1.pt",
        help="Path to second model checkpoint.",
    )
    parser.add_argument(
        "--model3",
        type=str,
        default="checkpoints/model_fold_2.pt",
        help="Path to third model checkpoint.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="predictions.csv",
        help="Output CSV path for batch predictions.",
    )
    parser.add_argument(
        "--n-top",
        type=int,
        default=50,
        help="Number of top compounds to return.",
    )
    parser.add_argument(
        "--host",
        type=str,
        default="0.0.0.0",
        help="Host for API serving mode.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port for API serving mode.",
    )
    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable auto reload in serving mode.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        choices=["auto", "cpu", "cuda"],
        help="Reserved for future device selection.",
    )
    return parser


def _run_predict(args: argparse.Namespace) -> int:
    if not args.input:
        print("Error: input H5AD file is required for batch prediction.", file=sys.stderr)
        return 1

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: input file '{args.input}' not found.", file=sys.stderr)
        return 1

    model_paths = [args.model1, args.model2, args.model3]
    missing = [path for path in model_paths if not Path(path).exists()]
    if missing:
        print(
            "Error: missing checkpoint files: " + ", ".join(missing),
            file=sys.stderr,
        )
        return 1

    adata = load_h5ad_file(str(input_path))
    model = DrugReflector(checkpoint_paths=model_paths)
    predictions = model.predict(adata, n_top=args.n_top)
    predictions.to_csv(args.output)

    print(f"Saved predictions to {args.output}")
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    if argv and argv[0] == "serve":
        from .api import app
        import uvicorn

        parser = _build_parser()
        args = parser.parse_args(argv[1:])
        uvicorn.run(app, host=args.host, port=args.port, reload=args.reload)
        return 0

    parser = _build_parser()
    args = parser.parse_args(argv)
    return _run_predict(args)


if __name__ == "__main__":
    raise SystemExit(main())
