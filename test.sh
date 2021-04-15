pylint src
pytest
mypy --disallow-untyped-defs src/
Rscript -e "lintr::lint_dir(path='src/', relative_path=TRUE)"
