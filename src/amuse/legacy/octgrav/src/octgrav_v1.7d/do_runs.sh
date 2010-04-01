mkdir ../compare/mw_$1
for theta in 0.8 0.7 0.6 0.5 0.4 0.3 0.2;
do
  echo $theta
  ./test_gravity ../MilkyWay/mw_$1.dat $theta  > ../compare/mw_$1/$theta.out 2> ../compare/mw_$1/$theta.log
done;
