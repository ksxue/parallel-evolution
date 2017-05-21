# This script puts the system to sleep until a certain job is done.
# It takes two arguments:
# 1: the job ID, name, or other identifier that is passed to -hold_jid
# 2: the length of time, in seconds, to sleep between checking whether the job is done

# Write a script wait.sh that will print "done" to a specified file.
echo 'echo "done" > $1' > wait.sh
chmod +x wait.sh

# Submit wait.sh to the queue.
# Put it on hold so that it finishes only after a job with the given name is done running.
qsub -hold_jid ${1} -cwd \
  -o wait.o -e wait.e \
  ./wait.sh ${1}.done

while [ ! -f ${1}.done ]
do
  sleep ${2}
done
rm -f ${1}.done
rm -f wait.sh
rm -f wait.o
rm -f wait.e
