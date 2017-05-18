#!/bin/bash

# Main repo hook file

REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../"
hookfile="$REPO_DIR/.git/hooks/pre-commit"

if [[ -f $hookfile ]]; then
    echo "'$hookfile' already exists - aborting" 1>&2
else
    echo '#!/bin/bash
    patch=$(git clang-format --diff)
    if [[ "$patch" =~ "no modified files to format" || "$patch" =~ "clang-format did not modify any files" ]]; then
       echo "" > /dev/null
    else
        echo ""
	echo "Your code has bad style." >&2
    	echo "Run 'git clang-format' to resolve the following issues:" >&2
	echo ""
    	echo "$patch"
    	exit 1
    fi
    ' > $hookfile
    chmod a+x $hookfile
fi
