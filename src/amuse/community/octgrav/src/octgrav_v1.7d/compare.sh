paste exact.out  tree > compare
../compare/process.rb < compare > dist_acc

# cat compare | awk '{print 1 - sqrt((($2-$4)/$2)^2)}'|sort -n| awk '{print 1-$1, (i++)/131072}' > dist_acc
# cat compare | awk '{print 1 - sqrt((($1-$3)/$1)^2)}'|sort -n| awk '{print 1-$1, (i++)/131072}' > dist_pot
