pylint src
pytest
mypy --disallow-untyped-defs src/
Rscript -e "lintr::lint_dir(path='.', relative_path=TRUE)"
