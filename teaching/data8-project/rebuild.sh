#!/bin/bash
#
# Regenerate the four notebooks from build_notebooks.py, then prove the two
# solutions notebooks still run top to bottom.
#
# Usage:  ./rebuild.sh          (from teaching/data8-project/)
#
# Only the solutions notebooks get executed. The student notebooks are full of
# `...` placeholders on purpose, so running those would fail by design.
#
# Executed output is thrown away (--stdout > /dev/null). That keeps the
# committed .ipynb files free of stored outputs, so their diffs stay readable.

# Run from this script's own directory, so it works from anywhere.
cd "$(dirname "$0")" || exit 1

failures=0

echo
echo "==> Regenerating notebooks from build_notebooks.py"
if python build_notebooks.py; then
    echo "    done"
else
    echo
    echo "!!! build_notebooks.py failed. The notebooks were not regenerated,"
    echo "!!! so nothing below would tell you anything. Fix the error above first."
    echo "!!! A stray indent inside an md(r\"\"\"...\"\"\") block is the usual cause."
    exit 1
fi

# Execute each solutions notebook. A clean run means every question's solution
# code still works against the current dataset and check_* helpers.
for track in datascience pandas; do
    notebook="solutions_${track}.ipynb"

    echo
    echo "==> Executing ${notebook}"

    log=$(mktemp)
    if jupyter nbconvert --to notebook --execute "$notebook" --stdout > /dev/null 2> "$log"; then
        echo "    clean"
    else
        echo
        echo "!!! ${notebook} failed to execute. Last 25 lines:"
        echo
        tail -n 25 "$log" | sed 's/^/    /'
        failures=$((failures + 1))
    fi
    rm -f "$log"
done

echo
if [ "$failures" -eq 0 ]; then
    echo "All good — 4 notebooks regenerated, both solutions run clean."
    echo "Open project_datascience.ipynb or project_pandas.ipynb to review."
else
    echo "${failures} of 2 solutions notebooks failed. See the output above."
fi
echo

exit "$failures"
