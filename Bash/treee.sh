tree -DFqs --timefmt "%F %T" --du \
| awk -v FS='[\]\[]' '{if ($2) {print "["$2"] "$1$3} else {print $0}}' 2>/dev/null \
| awk -v FS='──   '  '{if ($2) {print $1"── "$2}     else {print $0}}'

