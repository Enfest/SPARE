FOLDER="benches_unqomp/"

BENCHMARKS=(
 adder
 integercomparator
 wa
 mult
 wassaveqb
 mcry
 dj
 grover
 plr
)

# Vars
export PYTHONPATH=$(pwd)
export DOT_PATH="$(pwd)/dot_files/"
export LOG_PATH="$(pwd)/log_files/"
mkdir -p $(pwd)/"OUTPUT_JSON"/
mkdir -p $(pwd)/"OUTPUT_JSON"/"BENCH_RESULTS"
export json_folder=$(pwd)"/OUTPUT_JSON/BENCH_RESULTS/"

errors=()
RED='\033[0;31m'
NC='\033[0m' # No Color

basedir=$(pwd)
rewriter=$(pwd)/Spare/rewriter
export RIPL_HOME=$basedir

mkdir -p "$DOT_PATH"
mkdir -p $LOG_PATH
deactivate
cd ${basedir}/Unqomp
source /opt/conda/etc/profile.d/conda.sh
conda create --name unqomp --yes python=3.8
conda activate unqomp
pip install -e .
pip install cython=0.29.0

for b in ${!BENCHMARKS[@]}; do
	bench=${BENCHMARKS[$b]}
	json_path=$(json_folder)${bench}
	echo "Testing $bench..."
    python run_evaluation.py --file=json_path --bench=b
done

conda deactivate
cd ${basedir}
source venv/bin/activate

# python $basedir/src/bench_csv_aggregator.py $ $basedir/$benchout/suitesparse_check_sweep_reuse_${2}_${bench}.csv
echo -e "${RED}Failed tests:"
for i in ${!errors[@]}; do
    error=${errors[$i]} 
    echo -e "${RED}$error,"
done
echo -e "${NC}"
