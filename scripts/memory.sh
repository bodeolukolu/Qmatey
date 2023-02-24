#!/bin/bash
:> memory_usage.txt
while [[ ! -f "Analysis_Complete" ]]; do
  timestamp=$(date)
  usage=$(free -m | awk 'NR==2{print $3*100/$2}')
  if [[ "$usage" -gt 50 ]]; then
   echo $usage | awk -v pat="$timestamp" '{print $1"%\t"pat"\n"}' >> memory_usage.txt
  fi
  sleep 300
done
