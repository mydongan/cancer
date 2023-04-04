open(IN,$ARGV[0]);
my $line=<IN>;
chomp($line);
my @ll=split/\t/,$line;
my %hash;
my @LAG;
my @HAG;
my $gene;
my @info;
for(my $i=0;$i<=$#ll;$i++) {
	if ($ll[$i] =~ /^L.*tumor/ ) {push(@LAG,$i)}
	if ($ll[$i] =~ /^p.*tumor/ ) {push(@HAG,$i)}
	if ($ll[$i] =~ /Gene.refGene/ ) {$gene=$i}
	if ($ll[$i] =~ /^GO_ID|^GO_Term|^Pathway_ID|^Pathway_Des/) {push(@info,$i)}
}
print "Gene\tLAG.mut\tLAG.sum\tHAG.mut\tHAG.sum\tGO_ID\tGO_Term\tPathway_ID\tPathway_Des";
foreach my $i (@LAG) {print "\t$ll[$i]"}
foreach my $i (@HAG) {print "\t$ll[$i]"}
print "\n";
while(<IN>) {
	chomp;
	@ll=split/\t/;
	my $genes=$ll[$gene];
	my @ggs=split/,|;/,$genes;
	my %count;
	my @gg = grep { ++$count{ $_ } < 2; } @ggs;
	foreach my $g (@gg) {
		foreach my $i (@LAG) {
			if ($ll[$i] =~ /^\d.[1-9]/ ) {$hash{$g}{$i}++}
		}
		foreach my $i (@HAG) {
            if ($ll[$i] =~ /^\d.[1-9]/ ) {$hash{$g}{$i}++}
        }
		foreach my $i (@info) {
            $hash{$g}{info}.="\t$ll[$i]";
        }
	}
}
close IN;
foreach my $g (keys %hash) {
	print "$g";
	my @dats=split/\t/,$hash{$g}{info};
	my $data=join("\t",@dats[1..4]);
	my $ca=0;
	my $co=0;
	my $cas=@LAG;
	my $cos=@HAG;
	foreach my $i (@LAG) {
            if ( exists $hash{$g}{$i} ) {$ca++}
    }
	foreach my $i (@HAG) {
            if ( exists $hash{$g}{$i} ) {$co++}
    }
	print "\t$ca\t$cas\t$co\t$cos\t$data";
	foreach my $i (@LAG) {
            if ( exists $hash{$g}{$i} ) {print "\t$hash{$g}{$i}"} else {print "\t0"}
    }
	foreach my $i (@HAG) {
            if ( exists $hash{$g}{$i} ) {print "\t$hash{$g}{$i}"} else {print "\t0"}
    }
	print "\n";
}








