#/bin/sh

for f in $1/*.dat; do 
  f2=${f%.dat}
  python getWrithe.py $f $f2
done
