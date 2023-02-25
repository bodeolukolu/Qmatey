#!/bin/bash
:> memory_usage.txt
while [[ ! -f "Analysis_Complete" ]]; do
timestamp=$(date)
usage=$(free -m | awk 'NR==2{print $3*100/$2}')
echo $usage | awk -v pat="$timestamp" '{print $1"%\t"pat"\n"}'
sleep 300
done
