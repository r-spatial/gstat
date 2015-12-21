awk '{print "s/" $0 "/"}' qnames > tmp_a
sed 's/_/\\\\_/g' mapnames |  awk '{print "\\\\htmladdnormallink{" $0 "}"}' > tmp_b

awk '{print "{http:\/\/www.geog.uu.nl\/gstat\/manual\/gif\/" $0 ".gif}/g"}' mapnames > tmp_c

paste -d '' tmp_[abc]
rm tmp_[abc]
