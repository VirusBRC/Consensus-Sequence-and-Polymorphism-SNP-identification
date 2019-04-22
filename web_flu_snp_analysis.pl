#!/usr/bin/perl

#use strict;
use Cwd;
use Cwd 'chdir';
use Getopt::Std;
use lib ('/usr/lib64/R/library/RSPerl/perl/x86_64-linux-thread-multi');
#use R;


use vars qw($opt_r $opt_n $opt_a $opt_x $opt_t $opt_o $opt_g $opt_m $opt_c $opt_h $opt_d $opt_z $opt_e $opt_f $opt_s);
getopts('r:naxt:o:g:m:chdz:efs');

my $usage = "
     DESCRIPTION:
	This program do snp analysis. It will use muscle for protein seq and na sequence alignment if no pre-aligned na sequence available.


     OPTIONS:
	-r root directory for all data generated
	-n na seq only
	-a aa seq only 
	-x snp only 
	-t subtype 
	-o host
	-g segment
	-m protein (eg. PB1)
	-c run clustalw for na sequence alignment
	-h help
	-d debug flag
	-z min size
	-e reverse
	-f flu sequence
	-s skip uclust_muscle if input sequences already aligned	
";
		
if($opt_h) { die $usage; }

print "\t\tSTART: " . localtime(time) . "\n";

my $pwd = cwd();
if($opt_r) {$result_dir = $opt_r;}
&create_dir($result_dir);

my $main_log = $result_dir . "/main.log";
open(MAIN_LOG, ">$main_log") or die "cannot open log: $main_log to write $!\n";
print MAIN_LOG "\t\tSTART: " . localtime(time) . "\n\n";

my($na_only, $aa_only, $snp_only, $reverse, $isFlu, $skipMuscle) = (1,1,0,0,0,0);
if($opt_n){
	$aa_only = 0;
	print MAIN_LOG "Will only generate NA foma data\n";
}
if($opt_a){
	$na_only = 0;
	print MAIN_LOG "Will only generate AA foma data\n";
}
if($opt_x){
	$snp_only = 1;
	print MAIN_LOG "Will only generate NA SNP data\n";
}
if($opt_e){
	$reverse = 1;
	print MAIN_LOG "Will check reverse coding region\n";
}
if($opt_f){
	$isFlu = 1;
	print MAIN_LOG "Flu sequence\n";
}
if($opt_s){
	$skipMuscle = 1;
	print MAIN_LOG "Skip muscle\n";
}


my ($DEBUG);
if ($opt_d){$DEBUG = 1;}

my ($subtype, $host, $seqment, $prot, $min_size);
if($opt_t) {$subtype = $opt_t;}
if($opt_o) {$host = $opt_o;}
if($opt_g) {$seqment = $opt_g;}
if($opt_m) {$prot = $opt_m;}
if($opt_z) {$min_size = $opt_z;} else{$min_size = 250;}	

my $fulldataset = 0;
if ($subtype != undef){$fulldataset = 1;}

my $need_align = 0;
my $align_ok = 0;
if($na_only){
	my $na_file = $result_dir . "/" . "input.fasta";
	my $afa_file = $result_dir . "/" . "output.afa";
	my $aln_file = $result_dir . "/" . "output.aln";
	my $align_file = $result_dir . "/" . "alignment";
	my $gaf_file= $result_dir . "/" . "output.jpg";			
	if(!$snp_only){
		if($skipMuscle) {
			# "input.fasta" is aligned, just copy it over to the "output.afa" file
			`cp $na_file $afa_file`;
		} else {
			#my $cmd =  "muscle -in $na_file -fastaout $afa_file -clwout $aln_file";
			my $cmd =  "uclust_muscle -i $na_file -o $afa_file -a $align_file -c $aln_file";
			my $re = &run_command($cmd, 0);
			if($re != 0){
				print MAIN_LOG "\nMuscle failed on this group of $na_count flu DNA sequences from $subtype, $host, and $segment\n\n";
				print "\nMuscle failed on this group of $na_count flu DNA sequences from $subtype, $host, and $segment\n\n";
			}
		}
	}							
	
	#in case of input fasta is aligned, there will be no output.aln file and therefore no zip file will be generated
	if(-e $aln_file) {
		chdir $result_dir or die "cannot change to dir: $segment_dir\n";							
		my $aln_file = "output.aln";
		my $zip_file = "output.zip";
		my $cmd = "zip -r $zip_file $aln_file";
		&run_command($cmd, 1);
		chdir $pwd or die "cannot change to dir: $pwd\n";
	}
				
	my $log;
	
	
	($log, $gaf_file) = &generate_allele_freq($result_dir, $afa_file);			
	print MAIN_LOG "$log\n";		
		
	#($log) = &generate_snp($result_dir, $afa_file);			
	#print MAIN_LOG "$log\n";

	($log) = &generate_jpeg($gaf_file, $result_dir,$result_dir);
	print MAIN_LOG "$log\n\n";
		
	print MAIN_LOG "Finish snp analysis for na sequence\n\n";

} #end of na or snp loop

if($aa_only){ # start of aa loop
	my $aa_file = $result_dir . "/" ."input.fasta";
	my $afa_file  = $result_dir . "/" . "output.afa";
	my $aln_file = $result_dir . "/" . "output.aln";
	my $align_file = $result_dir . "/" . "alignment";
	if(!$snp_only){
		if($skipMuscle) {
			# "input.fasta" is aligned, just copy it over to the "output.afa" file
			`cp $aa_file $afa_file`;
		} else {						
			my $cmd = '';
			#$cmd = "muscle -in $aa_file -fastaout $afa_file -clwout $aln_file ";
			$cmd =  "uclust_muscle -i $aa_file -o $afa_file -a $align_file -c $aln_file";
			my $re = &run_command($cmd, 0);
			if($re != 0){
				print MAIN_LOG "\nMuscle failed on this group of $na_count flu protein sequences from $subtype, $host, $segment, and $prot\n\n";
				print "\nMuscle failed on this group of $na_count flu protein sequences from $subtype, $host, $segment, and $prot\n\n";
			}
		}
	}

	#in case of input fasta is aligned, there will be no output.aln file and therefore no zip file will be generated
	if(-e $aln_file){
		chdir $result_dir or die "cannot change to dir: $prot_dir\n";							
		my $aln_file = "output.aln";
		my $zip_file = "output.zip";
		my $cmd = "zip -r $zip_file $aln_file";
		&run_command($cmd, 1);
		chdir $pwd or die "cannot change to dir: $pwd\n";
	}

	my ($log) = &generate_protein_allele_freq($result_dir, $afa_file, $prot);
	print MAIN_LOG "$log\n\n";
		
	print MAIN_LOG "Finish snp analysis for aa sequence\n\n";	
	
} #end of aa loop

print MAIN_LOG "\t\tEND: " . localtime(time) . "\n";
close MAIN_LOG;

print "\t\tEND: " . localtime(time) . "\n";
exit(0);

sub generate_protein_allele_freq {
	my($prot_dir, $fasta, $prot) = @_;	

	my $log = '';
	my (%one2three) = qw(- Del A Ala C Cys D Asp E Glu F Phe G Gly H His I Ile K Lys L Leu M Met N Asn P Pro Q Gln R Arg S Ser T Thr V Val W Trp X Xaa Y Tyr);
	my $logBase = 2;
	my ($consFasta) = $prot_dir . "/cons.fasta";
	my ($fomaFile) = $prot_dir . "/foma.table";
	my (@conArray);
	my ($numSeq, $i) = (0, 0);
	my (%symbolHash);
	$/='>';
	open (FASTA, "$fasta") or die "can't open fasta_file $fasta\n";
	while (<FASTA>) {
    		chomp;
    		next unless (s/^([^\n]+)\n//);
    		$log .= "doing,$1\n";
    		s/\n//g; # remove \n from sequence
    		s/\r//g; # remove \r from sequence
    		
    		my (@seq) = split ('', $_);
    		my ($j, $nonGapStart, $nonGapEnd) = (0, 0, $#seq);
    		
    		$j++ while ($seq[$j] eq '-');
    		$nonGapStart = $j;
    		$j = $#seq;
    		$j-- while ($seq[$j] eq '-');
    		$nonGapEnd = $j;
    		for ($i=$nonGapStart;$i<=$nonGapEnd;$i++) {
        		
        		my ($symbol) = $seq[$i];
        		$symbolHash{$symbol} = 1;
        		if ($conArray[$i] && $conArray[$i]->{$symbol}) {
            			$conArray[$i]->{$symbol}++;
            			print "\n"
        		} else {
            			$conArray[$i]->{$symbol} = 1; #first time
        		}
    		}
    		$numSeq++;
    		
	}
#	return $log. "gap_start_in_all_seq\n" unless ($conArray[0]);
	my ($alignLen) = $#conArray + 1;
	$log .= "numSeq: $numSeq, alignLen: $alignLen\n";
	$i = $alignLen;
	my (@foma);
	my (@consensus);
	my (%unGap2GapMapping);
	my ($j, $k) = (0, 0);
	#@conArray=({'A'=>7,'T'=>1,'C'=>1,'G'=>1,'-'=>0},
	#           {'A'=>0,'T'=>8,'C'=>0,'G'=>0,'-'=>2},)
	while ($j <= $#conArray) {
    		my ($col) = $conArray[$j];
    		my ($conSymbol) = 'X';
    		my ($totalCount) = 0;
    		my (@symFreq);
    		foreach (values %$col) {$totalCount += $_;}
    		foreach my $symbol (keys %$col) {
        	my $symCount = $col->{$symbol};
        	push (@symFreq, $symCount/$totalCount);
        	$conSymbol = $symbol if ($symCount > $totalCount/2); 
    	}
	$foma[$j] = sprintf ("%.2f", &genProteinFOMA (\@symFreq, $logBase)) * 100;
    	$consensus[$j] = $conSymbol;
    	$conArray[$j]->{'totalSeq'} = $totalCount;
    	$conArray[$j]->{'consensus'} = $conSymbol;
    	$conArray[$j]->{'foma'} = $foma[$j];
    	if ($conSymbol ne '-') {
        	$unGap2GapMapping{$k} = $j;
        	$k++;
    	}
    		$j++;
	}
	$log .= "CS\t", join ("\t", @consensus), "\n";
	$log .= "FQ\t", join ("\t", @foma), "\n";
	my ($str) = join ('', @consensus);
	$str =~ s/-//g;
	open (CONS, ">$consFasta") or die "can't write cons fasta to current dir\n";
	print CONS ">consensus|$numSeq|$alignLen|", length $str, "\n$str\n";
	close CONS;

	open (FOMA, ">$fomaFile") or die "can't write to current dir $fomaFile\n";
	print FOMA "Position\tConsensus\tFOMA\tDetail\tNumberOfSequence\n";
	$i = 1;
	foreach my $col (@conArray) {
    		my (@detail);
    		while (my ($alphabet, $val) = each %$col) {
        		next unless (exists $symbolHash{$alphabet});
        		my ($three_letter) = $one2three{$alphabet};
        		push (@detail, "$three_letter=$val");
    		}
    		my ($cons1) = $col->{'consensus'}; #one letter aa code including del (-)
		$cons1 = $one2three{$cons1} unless ($cons1 eq '-');
            	my $position = $i;
                if ($cons1 eq '-'){
                        $position = 'N/A';
                        $i--;
                } 
		if (!$fulldataset){
	    		print FOMA "$position\t$cons1\t", $col->{'foma'}, "\t",
	            	join (',', sort @detail), "\t", $col->{'totalSeq'}, "\n";
		}else{	
			print FOMA "$subtype\t$host\t$segment\t$prot\t$position\t$cons1\t", $col->{'foma'}, "\t",
	            	join (',', sort @detail), "\t", $col->{'totalSeq'}, "\n";
		}
    		$i++;
	}
	close FOMA;
}

sub genProteinFOMA {
    my ($symFreq, $logBase) = @_;
    my ($numSym) = scalar @$symFreq;
    my ($foma, $i, $totalFreq) = (0, 0, 0);
    for($i=0;$i<$numSym;$i++) {
        my ($pi) = $symFreq->[$i];
        $totalFreq += $pi;
        $foma += $pi * (log($pi)/log($logBase));
    }
    die "cumulative frequency > 1 $totalFreq\n" if (sprintf ("%.6f", $totalFreq) > 1);
    #return log(4)/log(2) + $foma;
    return 0 - $foma;
}

sub generate_snp{
	my ($segment_dir, $fasta) = @_;

	my $snp_file = $segment_dir . "/snp.out";
	my $log = '';
	my (@conArray);
	my ($numSeq, $i) = (0, 0);
	my (%acc2aln);
	$/='>';
	my $out = 1;
	open (FASTA, "$fasta") or die "can't open fasta_file $fasta\n";
	while (<FASTA>) {
    		chomp;
    		next unless (s/^([^\n]+)\n//);
    		my ($header) = $1;
		$log .= "doing,$header\n";
    		#return  $log . "problem_with_header $header\n" unless ($header =~ /^([^\|]+)\|/);
	    	my ($acc) = $1;
		
#		if($acc eq "DQ021659"){print "$_\n";print length($_);print "\n";}
#		if($acc eq "AY699987"){print "$_\n";print length($_);print "\n";}
#		print "$acc  _LEN_" . length($_);
#		print "\n";  		

		s/\n//g; # remove \n from sequence
		s/\r//g; # remove \r from sequence
    		my (@seq) = split ('', $_);
    		$acc2aln{$acc} = \@seq;
    		my ($j, $nonGapStart, $nonGapEnd) = (0, 0, $#seq);
    		$j++ while ($seq[$j] eq '-');
    		$nonGapStart = $j;
    		$j = $#seq;
		$j-- while ($seq[$j] eq '-');
    		$nonGapEnd = $j;
    		for ($i=$nonGapStart;$i<=$nonGapEnd;$i++) {
			$out = 0;
        		my ($symbol) = $seq[$i];
        		if ($conArray[$i] && $conArray[$i]->{$symbol}) {
            			$conArray[$i]->{$symbol}++;
        		} else {
            			$conArray[$i]->{$symbol} = 1; #first time
        		}
    		}
    		$numSeq++;
	}
	close FASTA;
	return $log . "gap_start_in_all_seq\n" if ($out);
	my ($alignLen) = $#conArray + 1;
	$log .= "numSeq: $numSeq, alignLen: $alignLen\n";
	my (@consensus);
	my ($j) = 0;
	#@conArray=({'A'=>7,'T'=>1,'C'=>1,'G'=>1,'-'=>0},
	#           {'A'=>0,'T'=>8,'C'=>0,'G'=>0,'-'=>2},)
	while ($j <= $#conArray) {
    		my ($col) = $conArray[$j];
    		my ($conSymbol) = 'N';
    		my ($totalCount) = 0;
    		my (@symFreq);
    		foreach (values %$col) {$totalCount += $_;}
    		foreach my $symbol (keys %$col) {
        		my $symCount = $col->{$symbol};
        		push (@symFreq, $symCount/$totalCount);
        		$conSymbol = $symbol if ($symCount > $totalCount/2); 
    		}
    		$consensus[$j] = $conSymbol;
    		$conArray[$j]->{'totalSeq'} = $totalCount;
    		$conArray[$j]->{'consensus'} = $conSymbol;
    		$j++;
	}

	#add code to remove start trash seq
	my $array_len = @conArray;
	my $start_pos = 0;
	my $end_pos = $array_len -1;
	if ($isFlu){
	for (my $i=0; $i<$array_len; $i++){
		my $str = $conArray[$i]->{'consensus'}.$conArray[$i+1]->{'consensus'}.$conArray[$i+2]->{'consensus'};
		if($str eq "AGC"){$start_pos = $i; last;}
		if($conArray[$i]->{'totalSeq'} > 10){ last;}
		if($i > 40) {last;}
	}

	for (my $i=$array_len-1; $i>0; $i--){
		my $str = $conArray[$i-4]->{'consensus'}.$conArray[$i-3]->{'consensus'}.$conArray[$i-2]->{'consensus'}.$conArray[$i-1]->{'consensus'}.$conArray[$i]->{'consensus'};
		if($str eq "CTACT"){$end_pos = $i;last;}
		if($conArray[$i]->{'totalSeq'} > 10){ last;}
		if($i < $array_len-40) {last;}
	}
	}
	$log .= "LEN: $array_len STAT: $start_pos   END: $end_pos\n";


	#k  :    1234   5678  90              index of seq (w/o gaps)
	#j  : 0123456789012345678901          index of aln and cons
	#aln: ---ATCG---TTGG--AC---
	#con: CCTAACGTT-CA---TACCT--
	open(SP, ">$snp_file") or die "cannot open snp_file: $snp_file to write $!\n";
	while (my ($acc, $aln) = each %acc2aln) {
		#this need rework
#		print "$subtype"."_"."$host"."_"."$segment"."problem_with_cons_aln $acc" ."_". (scalar @consensus) ."_".  (scalar @$aln)."\n" unless ((scalar @consensus) == (scalar @$aln));
	 	return $log . "problem_with_cons_aln $acc\n" unless ((scalar @consensus) == (scalar @$aln));
	 	#    if ($acc eq 'L06579') {
		#    print STDERR  join('', @consensus), "\n";
		#    print STDERR  join('', @$aln), "\n";}
    		my ($n, $m, $k) = (0, $#consensus, 1);
    		#comparison starts after leading -
    		$n++ while ($aln->[$n] eq '-');
    		my $l=1;
    		if($n < $start_pos) {$n = $start_pos; } else {$l = $n+1;}
    		#comparison ends before trailing -
    		$m-- while ($aln->[$m] eq '-');
    		if($m > $end_pos) {$m = $end_pos;}
    		for (my $i=0; $i<$n; $i++){
			if($aln->[$i] ne '-'){
				$k++;
			}	
    		}

    		for ($j = $n; $j <= $m; $j++) {
        		my ($symCon, $symAcc) = ($consensus[$j], $aln->[$j]);
        		if ($symAcc eq '-' && $symCon ne '-') { #deletion
				if (!$fulldataset){
            				print SP "$acc\t", $k-1, "\t", $l, "\tdeletion\t$symCon\n";
				} else {
					print SP "$subtype\t$host\t$segment\t$acc\t", $k-1, "\t", $l, "\tdeletion\t$symCon\n";
				}
				
		        } elsif ($symCon eq '-' && $symAcc ne '-') { #insertion
				if (!$fulldataset){
            				print SP "$acc\t$k\t", $l, "\tinsertion\t$symAcc\n";
				} else {
					print SP "$subtype\t$host\t$segment\t$acc\t$k\t", $l, "\tinsertion\t$symAcc\n";
				}
				
			      	$k++;
        		} elsif ($symCon eq '-' && $symAcc eq '-') { #match of gaps
		            	#do nothing
        		} else { #either identical bases or mismatch
				if (!$fulldataset){
            				print SP "$acc\t$k\t", $l, "\tmismatch\t$symCon" . '->' . "$symAcc\n" if ($symCon ne $symAcc); #mismatch
				} else {
					print SNP "$subtype\t$host\t$segment\t$acc\t$k\t", $l, "\tmismatch\t$symCon" . '->' . "$symAcc\n" if ($symCon ne $symAcc); #mismatch
				}
				
            			$k++;
        		}
			$l++;
    		}
	}	
	close SP;
	return ($log);
}

sub generate_jpeg{
	my($gaf_file, $segment_dir, $na_jpg_dir) = @_;

	$name = 'Rplot'; # unless ($name);
	$type = 'jpeg'; # unless ($type);

#	&R::startR("");
#	&R::startR("--silent");
	return unless (-s "$gaf_file");
	
	my $pwd = cwd();
    	chdir $segment_dir;
    	my ($cmd) = "grep '^FQ' gaf.log";
    	my $reStr = `$cmd`;
    	$reStr =~ s/FQ//;
    	$reStr =~ s/\n//;
	return "grep FQ gaf.log failed $dir\n" unless ($reStr =~ /^(\d|\t)+$/);
#    die "grep FQ gaf.log failed $dir\n" unless ($reStr =~ /^(\d|\t)+$/);
	my (@fqs) = split (/\t/, $reStr);
    	my ($len) = scalar @fqs;
    	my (@x) = 1..$len;
    	if ($type eq 'jpeg') {
#       	&R::jpeg ("$name.$type"); 
    	} else {
#        	&R::pdf ("$name.$type"); 
    	}
#   	&R::callWithNames("plot", {"x"=>\@x, "y"=>\@fqs, "xlab"=>"position", "ylab"=>"score", "type"=>"l"}); 
#    	&R::call("dev.off"); 
	
	chdir $pwd;
	return;
}


sub generate_allele_freq {
	my ($data_dir, $fasta) = @_;
	my $log = '';
	my ($consFasta) = $data_dir . "/cons.fasta";
	my ($logBase) = 2;
	my ($getCDScmd) = '/usr/local/bin/getorf -outseq stdout -find 3';
	if (!$reverse) {$getCDScmd = $getCDScmd." -noreverse";}
	$getCDScmd = $getCDScmd." 2> /dev/null -minsize ".$min_size;
	my ($fomaFile) = $data_dir . "/foma.table";
	my ($gaf_file) = $data_dir . "/gaf.log";
	open(GF, ">$gaf_file") or die "cannot open file: $gaf_file to write $!\n";
	my (@conArray);
	my ($numSeq, $i) = (0, 0);
	my (%symbolHash);
	$/='>';
	my $out = 1;
	open (FASTA, "$fasta") or die "can't open fasta_file $fasta\n";
	while (<FASTA>) {
    		chomp;
    		next unless (s/^([^\n]+)\n//);
    		print GF "doing,$1\n";
    		s/\n//g; # remove \n from sequence
    		s/\r//g; # remove \r from sequence
		#change lowercase nucleotides to uppercase
		$_ =~ tr/atgc/ATGC/;
    		my (@seq) = split ('', $_);
    		my ($j, $nonGapStart, $nonGapEnd) = (0, 0, $#seq);
    		$j++ while ($seq[$j] eq '-');
    		$nonGapStart = $j;
    		$j = $#seq;
    		$j-- while ($seq[$j] eq '-');
    		$nonGapEnd = $j;
    		for ($i=$nonGapStart;$i<=$nonGapEnd;$i++) {
			$out = 0;
        		my ($symbol) = $seq[$i];
			if($symbol eq 'N'){
				next;
			}
        		$symbolHash{$symbol} = 1;
        		if ($conArray[$i] && $conArray[$i]->{$symbol}) {
            			$conArray[$i]->{$symbol}++;
        		} else {
            			$conArray[$i]->{$symbol} = 1; #first time
        		}
    		}
    		$numSeq++;
	}
	close FASTA;
# this may be need to be added
	return $log if ($out);
#	return $log unless ($conArray[0]);
#	die "gap_start_in_all_seq\n" unless ($conArray[0]);

	my ($alignLen) = $#conArray + 1;
	print GF "numSeq: $numSeq, alignLen: $alignLen\n";
	my (%printHash);
	$i = $alignLen;
	# initialize %printHash=('A'=>[0,0,0,...,0], 'T'=>[0,0,0,...,0],
	#                       'C'=>[0,0,...,0], 'G'=>[0,0,0,...,0], '-'=>	[0,0,0,...,0])
	foreach my $alphabet (keys %symbolHash) {
    		while ($i >= 0) {
        		$i--;
        		$printHash{$alphabet}->[$i] = 0; 
    		}
    		$i=$alignLen;
	}
	my (@foma);
	my (@consensus);
	my (%unGap2GapMapping);
	my ($j, $k) = (0, 0);
	#@conArray=({'A'=>7,'T'=>1,'C'=>1,'G'=>1,'-'=>0},
	#           {'A'=>0,'T'=>8,'C'=>0,'G'=>0,'-'=>2},)
	while ($j <= $#conArray) {
    		my ($col) = $conArray[$j];
    		my ($conSymbol) = 'N';
    		my ($totalCount) = 0;
    		my (@symFreq);
    		foreach (values %$col) {$totalCount += $_;}
    			foreach my $symbol (keys %$col) {
        			my $symCount = $col->{$symbol};
        			$printHash{$symbol}->[$j] = $symCount;
        			push (@symFreq, $symCount/$totalCount);
        			$conSymbol = $symbol if ($symCount > $totalCount/2); 
    			}
    			$foma[$j] = sprintf ("%.2f", &genFOMA (\@symFreq, $logBase)) * 100;
	    		$consensus[$j] = $conSymbol;
    			$conArray[$j]->{'totalSeq'} = $totalCount;
    			$conArray[$j]->{'consensus'} = $conSymbol;
    			$conArray[$j]->{'foma'} = $foma[$j];
    			if ($conSymbol ne '-') {
        			$unGap2GapMapping{$k} = $j;
        			$k++;
    			}
    			$j++;
		}
		foreach my $alphabet (sort keys %printHash) {
    			print GF "$alphabet\t", join ("\t", @{$printHash{$alphabet}}), "\n";
		}
		print GF "CS\t", join ("\t", @consensus), "\n";
		print GF "FQ\t", join ("\t", @foma), "\n";
		my ($str) = join ('', @consensus);
		$str =~ s/-//g;
		open (CONS, ">$consFasta") or die "can't write cons fasta to current dir\n";
		print CONS ">consensus|$numSeq|$alignLen|", length $str, "\n$str\n";
		close CONS;
		my ($coords) = &getCDScoords ($getCDScmd, $consFasta);
#		return ($log, $gaf_file) unless (scalar @$coords);
		goto LABEL unless (scalar @$coords);
		my (@unGappedCDScoords) = @$coords;

#check the leftmost begin coord is ATG and rightmost end coord + 3 is STOP
		my (@sortBegin) = sort {$a->{'begin'} <=> $b->{'begin'}} @unGappedCDScoords;
		my (@sortEnd) = sort {$a->{'end'} <=> $b->{'end'}} @unGappedCDScoords;
		my ($leftmost, $rightmost) = ($sortBegin[0]->{'begin'}, $sortEnd[-1]->{'end'});
		my ($firstCodon) = substr ($str, $leftmost - 1, 3);
#		print STDERR "no_ATG_coding_start\n" unless ($firstCodon eq 'ATG');
		my ($lastCodon) = substr ($str, $rightmost, 3);
#print STDERR "no_STOP_codon\n" unless ($lastCodon =~ /^T(AA|AG|GA)$/);
		unless ($lastCodon =~ /^T(AA|AG|GA)$/) { #no STOP and 3'UTR
#			print STDERR "no_STOP_codon\n";
		$rightmost--; #rightmost was the coord of STOP, now it's the end of seq
	}
	LABEL:
	my (@maskedArray);
	for ($i=0; $i<=$#consensus; $i++) { #initialize to non-coding
    		$maskedArray[$i] = 0;
	}
	#set coding region to 1
	my ($gappedCoordBegin) = $unGap2GapMapping{$leftmost};
	my ($gappedCoordEnd) = $unGap2GapMapping{$rightmost};
	my ($n) = 0;
	for ($n = $gappedCoordBegin - 1; $n < $gappedCoordEnd; $n++) {
		foreach my $coords (@unGappedCDScoords) {
    			if (($n > $coords->{'begin'} && $n < $coords->{'end'}) ||($n > $coords->{'end'} && $n < $coords->{'begin'}) ){
				$maskedArray[$n] = 1;
			}
		}    
	}

	my $array_len = @conArray;
	my $start_pos = 0;
	my $end_pos = $array_len -1;
	if ($isFlu){
	for (my $i=0; $i<$array_len; $i++){
		my $str = $conArray[$i]->{'consensus'}.$conArray[$i+1]->{'consensus'}.$conArray[$i+2]->{'consensus'};
		if($str eq "AGC"){$start_pos = $i; last;}
		if($conArray[$i]->{'totalSeq'} > 10){ last;}
		if($i > 40) {last;}
	}

	for (my $i=$array_len-1; $i>0; $i--){
		my $str = $conArray[$i-4]->{'consensus'}.$conArray[$i-3]->{'consensus'}.$conArray[$i-2]->{'consensus'}.$conArray[$i-1]->{'consensus'}.	$conArray[$i]->{'consensus'};
		if($str eq "CTACT"){$end_pos = $i;last;}
		if($conArray[$i]->{'totalSeq'} > 10){ last;}
		if($i < $array_len-40) {last;}
	}
	}
	#print "LEN: $array_len STAT: $start_pos   END: $end_pos\n";

	open (FOMA, ">$fomaFile") or die "can't write to current dir $fomaFile\n";
	print FOMA "Position\tCoding\tFOMA\tConsensus\tA\tT\tG\tC\tDeletion\tNumberOfSequence\n";
	$i = 1;
	for($j=$start_pos; $j<=$end_pos; $j++){
	#foreach my $col (@conArray) {
    		my $col = $conArray[$j];
    		foreach my $alphabet (keys %symbolHash) {
        		$col->{$alphabet} = 0 unless (exists $col->{$alphabet});
    		}
    # when there is no gap in the alignment, '-' won't showup in %symbolHash
    # and $col->{'-'} undefined. need to set it to 0.
    		unless (exists $col->{'-'}) {$col->{'-'} = 0;}
    		my ($coding) = 'yes';
		my $position = $i;
    		if ($col->{'consensus'} eq '-') {
        		$coding = 'N/A';
			$position = 'N/A';
			$i--;
    		} elsif ($maskedArray[$i - 1] == 0) { #$i starts at 1
        		$coding = 'no';
    		}
		if (!$fulldataset){
		    	print FOMA "$position\t$coding\t", $col->{'foma'}, "\t", $col->{'consensus'},
		        "\t", $col->{'A'}, "\t", $col->{'T'}, "\t", $col->{'G'},
		        "\t", $col->{'C'}, "\t", $col->{'-'}, "\t", $col->{'totalSeq'}, "\n";
		} else {
			print FOMA "$subtype\t$host\t$segment\t$position\t$coding\t", $col->{'foma'}, "\t", $col->{'consensus'},
		        "\t", $col->{'A'}, "\t", $col->{'T'}, "\t", $col->{'G'},
		        "\t", $col->{'C'}, "\t", $col->{'-'}, "\t", $col->{'totalSeq'}, "\n";
		}
    		$i++;
	}
	close FOMA;
	close GF;

	return ($log, $gaf_file);
}

sub genFOMA {
    my ($symFreq, $logBase) = @_;
    my ($numSym) = scalar @$symFreq;
    my ($foma, $i, $totalFreq) = (0, 0, 0);
    for($i=0;$i<$numSym;$i++) {
        my ($pi) = $symFreq->[$i];
        $totalFreq += $pi;
        $foma += $pi * (log($pi)/log($logBase));
    }
    die "cumulative frequency > 1 $totalFreq\n" if (sprintf ("%.6f", $totalFreq) > 1);
    #return log(4)/log(2) + $foma;
    return 0 - $foma;
}

sub getCDScoords {
    my ($getCDScmd, $consFile) = @_;
    my (@CDScoords);
    my ($cmd) = "$getCDScmd $consFile | grep '^>'";
   
    my (@reStr) = `$cmd`;
    foreach (@reStr) {
        next unless (/\[(\d+) - (\d+)\]/);
#        print "CDScoords ", $_;
        push (@CDScoords, {'begin' => $1, 'end' => $2});
    }
#    die "getCDScmd failed\n" unless (scalar @CDScoords);
    return \@CDScoords;
}

sub run_command{
	my($cmd, $fatal) = @_;

	if($cmd){
		my $re = system($cmd);
		
		if($fatal && ($re != 0)){
			print MAIN_LOG "Fatal error: cannot excute cmd: $cmd\n";
			print "\nFatal error: cannot excute cmd: $cmd\n";
			exit(-1);	
		} elsif (!$fatal && ($re != 0)){
			print MAIN_LOG "Waring: cannot excute cmd: $cmd\n";
			print "\nWarning: cannot excute cmd: $cmd\n";
			return (-1);
		}
	}
	return 0;
}

sub create_dir{
	my ($dir) = @_;
		
	if(!(-d $dir)) {	
		mkdir($dir) || die "Failed making $dir<br>\n";
		chmod(0777,$dir) || print "Failed chmod $dir<br>\n";
	}
}







