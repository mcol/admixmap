#!/bin/bash
# 
# This scripts runs other scripts: submits the tasks, waits for them to
# finish, processes the output, uploads results to web, sends an e-mail
# notification and serves coffee.

CONFIG_FILE="analyses-20"
# 
# bash manage-analyses.sh \
# 	--action setup \
# 	--config "${CONFIG_FILE}" \
# 	--submit \
# 	--data /ichec/home/users/maciej/shared-genepi/maciej/chr22-tuned-2/hapmap
# 
# sleep 10

echo "Waiting for the tasks to finish"
while [ ! -z "$(qstat | grep maciej)" ]
do
	sleep 60
done

bash manage-analyses.sh \
	--action results \
	--config "${CONFIG_FILE}"

rsync -az html/* wahwah@suita.chopin.edu.pl:/home/student/wahwah/html/ucd

(echo "Howdy Maciej, your analyses have finished."; \
echo "Please visit http://suita.chopin.edu.pl/~wahwah/ucd/ for results.") | mail -s "Analyses finished" maciej.blizinski@ucd.ie -R maciej.blizinski@ucd.ie
