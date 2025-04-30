FOLDER="benches_length_simplified/"

BENCHMARKS=(
 length_simplified_orig1
 length_simplified1
 length_simplified_only_cn1
 length_simplified_only_cf1
 length_simplified5
 length_simplified_only_cn5
 length_simplified_only_cf5
 length_simplified_orig5
 length_simplified9
 length_simplified_orig9
 length_simplified_only_cf9
 length_simplified_only_cn9
)

# Vars
export PYTHONPATH=$(pwd)
export DOT_PATH="$(pwd)/dot_files/"
export LOG_PATH="$(pwd)/log_files/"

errors=()
RED='\033[0;31m'
NC='\033[0m' # No Color

basedir=$(pwd)
rewriter=$(pwd)/Spare/rewriter
export RIPL_HOME=$basedir


mkdir -p "$DOT_PATH"
mkdir -p $LOG_PATH

for b in ${!BENCHMARKS[@]}; do
	bench=${BENCHMARKS[$b]}
	mkdir -p $DOT_PATH/$FOLDER/${bench}
	mkdir -p $LOG_PATH/$FOLDER/${bench}
	mkdir -p $(pwd)/"OUTPUT_JSON"/
	mkdir -p $(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS"
	mkdir -p $(pwd)/"OUTPUT_JSON"/"BENCH_RESULTS"
	mkdir -p $(pwd)/"temp_files"/
	dotpath=$DOT_PATH/$FOLDER/${bench}/
	rewritelogpath=$LOG_PATH/$FOLDER/${bench}/"rewritelog5.txt"
	manuallogpath=$LOG_PATH/$FOLDER/${bench}/"benchlog.txt"
	tool_json_path=$(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS"
	bench_json_path=$(pwd)/"OUTPUT_JSON"/"BENCH_RESULTS"
	echo "Testing $bench..."
	echo $rewriter/compile_rewrite_graphs.py
	python $rewriter/compile_rewrite_graphs.py --test-name=$dotpath${bench}_graph --compile-mode="optimized_qubit"\
		--runs=$1 --json-name=$tool_json_path/output_$bench.json > $rewritelogpath
	# python $rewriter/compile_rewrite_graphs.py --test-name=$dotpath${bench}_graph \
	# 	--compile-mode="bench_qubit" --runs=$1 --json-name=$bench_json_path/output_$bench.json > $manuallogpath
	rm -r $(pwd)/"temp_files"/
done

# python $basedir/src/bench_csv_aggregator.py $ $basedir/$benchout/suitesparse_check_sweep_reuse_${2}_${bench}.csv
echo -e "${RED}Failed tests:"
for i in ${!errors[@]}; do
    error=${errors[$i]} 
    echo -e "${RED}$error,"
done
echo -e "${NC}"
