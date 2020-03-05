source /software/intel/advisor/advixe-vars.sh intel64

opts="--project-dir=./advi --search-dir all:=../src/ -- ./cato.x input.ini"

# Collect the data
advixe-cl --collect=survey ${opts}

# Save the results
advixe-cl --report=survey --project-dir=./advi --format=text --report-output=./out/survey.txt
