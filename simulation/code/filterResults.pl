#!/usr/env/ perl
use warnings;
use strict;
use FileHandle;
use Data::Dumper;

#input parameter file index
my ($paramIndex) = @ARGV;
my $resultsDir = "resultsDir";

my %paramResults;

my $paramFile = "$resultsDir\/param_ALLDRUG_$paramIndex.txt";
my $PARAM = FileHandle->new("$paramFile") or die("$! $paramFile");

my $stoptFile = "$resultsDir\/stopt_ALLDRUG_$paramIndex.txt";
my $STOPT = FileHandle->new("$stoptFile") or die("$! $stoptFile");

my $dosageFile = "$resultsDir\/dosage_ALLDRUG_$paramIndex.txt";
my $DOSAGE = FileHandle->new("$dosageFile") or die("$! $dosageFile");

#param file
#(1)parameter configuration index (one entry).
#(2)initial population compositions x0 (four entries, x0(S), x0(R1), x0(R2), x0(R12)).
#(3)growth rate g0 (one entry).
#(4)drug sensitivity matrix Sa (eight entries, Sa(S,drug1), Sa(S,drug2), Sa(R1,drug1), Sa(R1,drug2), ...).
#(5)transition rate matrix T (sixteen entries, T(S->S), T(S->R1), T(S->R2), T(S->R12), T(R1->S), T(R1->R1), T(R1->R2), T(R1->R12)
foreach my $paramLine (<$PARAM>){
	$paramLine =~s/\R//;
	my @paramValues = split(/ /,$paramLine);
	my $paramID = $paramValues[0];
	my $g0 = $paramValues[5];
	my $S_D1 = $paramValues[6];
	my $R12_D1 = $paramValues[12];
	my $R12_D2 = $paramValues[13];
	$paramResults{$paramID}{"g0"} = $g0;
	$paramResults{$paramID}{"S_D1"} = $S_D1;
	$paramResults{$paramID}{"R12_D1"} = $R12_D1;
	$paramResults{$paramID}{"R12_D2"} = $R12_D2;	
	
	if ($S_D1 > 0.2*$g0){
		$paramResults{$paramID}{"paramFilter"} = 1;
	}else{
		$paramResults{$paramID}{"paramFilter"} = 0;
	}

	if (($R12_D1 < $g0) && ($R12_D2 < $g0)){
		$paramResults{$paramID}{"resistantFilter"} = 1;
	}else{
		$paramResults{$paramID}{"resistantFilter"} = 0;
	}
		
}

#stopt file
#strategy 0: strategy 0 in the PNAS paper.
#strategy 1: strategy 1 in the PNAS paper.
#strategy 2: strategy 2.1 in the PNAS paper.
#strategy 3: strategy 2.2 in the PNAS paper.
#strategy 4: strategy 3 in the PNAS paper.
#strategy 5: dynamic programming extension of strategy 1 in the Biology Direct paper.
#strategy 6: dynamic programming extension of strategy 2.1 in the Biology Direct paper.
#strategy 7: dynamic programming extension of strategy 2.2 in the Biology Direct paper.
#strategy 8: dynamic programming extension of strategy 3 in the Biology Direct paper.
#strategy 9: global dynamic programming in the Biology Direct paper.
foreach my $stoptLine (<$STOPT>){
	$stoptLine =~s/\R//;
	my @stoptValues = split(/ /,$stoptLine);
	my $paramID = $stoptValues[0];
	my $stopt_strgy0 = $stoptValues[1];
	my $stopt_strgy22 = $stoptValues[2];
	my $stopt_strgy22trial = $stoptValues[3];
	$paramResults{$paramID}{"stopt_Strgy0"} = $stopt_strgy0;
	$paramResults{$paramID}{"stopt_Strgy2.2"} = $stopt_strgy22;
	$paramResults{$paramID}{"stopt_Strgy2.2trial"} = $stopt_strgy22trial;		
	if ($stopt_strgy0 <= 1460){
		$paramResults{$paramID}{"stoptFilter"} = 1;
	}else{
		$paramResults{$paramID}{"stoptFilter"} = 0;
	}
}

#dosage file
#(1)parameter configuration index.
#(2)strategy index (0-9).
#(3)(drug1 dosage,drug2 dosage) at t=0.
#(4)(drug1 dosage,drug2 dosage) at t=45.
#...
#(42)(drug1 dosage,drug2 dosage) at t=1755.
foreach my $dosageLine (<$DOSAGE>){
	$dosageLine =~ s/\R//;
	my @dosageValues = split(/ /,$dosageLine);
	my $strategy = $dosageValues[1];
	if ($strategy == 0){
		my $paramID = $dosageValues[0];
		my @dosage_strgy0 = @dosageValues[2,3,4,5];
		$paramResults{$paramID}{"dosage_strgy0"} = [ @dosage_strgy0 ];
	}
	if ($strategy == 1){
		my $paramID = $dosageValues[0];
		my @dosage_strgy22 = @dosageValues[2,3,4,5];
		$paramResults{$paramID}{"dosage_strgy2.2"} = [ @dosage_strgy22 ];
	}
	if ($strategy == 2){
		my $paramID = $dosageValues[0];
		my @dosage_strgy22trial = @dosageValues[2,3,4,5];
		$paramResults{$paramID}{"dosage_strgy2.2trial"} = [ @dosage_strgy22trial ];
	}
}



foreach my $paramID (keys %paramResults){
	if (
		($paramResults{$paramID}{"paramFilter"} == 1) && 
		($paramResults{$paramID}{"stoptFilter"} == 1) &&
		($paramResults{$paramID}{"resistantFilter"} == 1)
		)
	{
		my @moves_strgy0 = @{ $paramResults{$paramID}{"dosage_strgy0"} };
		my @moves_strgy22 = @{ $paramResults{$paramID}{"dosage_strgy2.2"} };

		if (($moves_strgy0[0] == $moves_strgy22[0]) &&
			($moves_strgy0[1] == $moves_strgy22[1]))
		{
			if (($moves_strgy0[2] == $moves_strgy22[2]) &&
				($moves_strgy0[3] == $moves_strgy22[3]))
			{
				$paramResults{$paramID}{"dosageCompare"} = "bothSame";
			}else{
				$paramResults{$paramID}{"dosageCompare"} = "firstSame";
			}
		}else{
			if (($moves_strgy0[2] == $moves_strgy22[2]) &&
				($moves_strgy0[3] == $moves_strgy22[3]))
			{
				$paramResults{$paramID}{"dosageCompare"} = "secondSame";
			}else{
				$paramResults{$paramID}{"dosageCompare"} = "bothDiff";
			}
		}
		my $compare = $paramResults{$paramID}{"dosageCompare"};
		my $stopt_strgy0 = $paramResults{$paramID}{"stopt_Strgy0"};
		my $stopt_strgy22 = $paramResults{$paramID}{"stopt_Strgy2.2"};
		my $stopt_strgy22trial = $paramResults{$paramID}{"stopt_Strgy2.2trial"};
		print "$paramID\t$compare\tStrgy0\t$stopt_strgy0\n";
		print "$paramID\t$compare\tStrgy2.2\t$stopt_strgy22\n";
		print "$paramID\t$compare\tStrgy2.2trial\t$stopt_strgy22trial\n";
		#print Dumper($paramResults{$paramID});
		#print Dumper($paramID);
		#print Dumper(@moves_strgy0,@moves_strgy22);
	}
	#print Dumper($paramResults{$paramID});
}
