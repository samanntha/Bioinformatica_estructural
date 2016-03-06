#!/usr/bin/perl -w

# Programa: Predicción de promotores
# Licenciatura en Ciencias Genómicas
# Integrantes del equipo: Jessica Samantha Cruz Ruiz, Lorena Elizabeth Fajardo Brigido

use strict;
use warnings;

# global variables
my $T           = 37; # temperature(C)
my $windowL     = 15;  # window length, http://www.biomedcentral.com/1471-2105/6/1
my %NNparams    = ( 
	# SantaLucia J (1998) PNAS 95(4): 1460-1465.
	# [NaCl] 1M, 37C & pH=7 
	# H(enthalpy): kcal/mol	, S(entropy): cal/k�mol
	# stacking dinucleotides
	'AA/TT' , {'H',-7.9, 'S',-22.2},
	'AT/TA' , {'H',-7.2, 'S',-20.4},
	'TA/AT' , {'H',-7.2, 'S',-21.3},
	'CA/GT' , {'H',-8.5, 'S',-22.7},
	'GT/CA' , {'H',-8.4, 'S',-22.4},
	'CT/GA' , {'H',-7.8, 'S',-21.0},
	'GA/CT' , {'H',-8.2, 'S',-22.2},
	'CG/GC' , {'H',-10.6,'S',-27.2},
	'GC/CG' , {'H',-9.8, 'S',-24.4},
	'GG/CC' , {'H',-8.0, 'S',-19.9},
	# initiation costs
	'G'     , {'H', 0.1, 'S',-2.8 },
	'A'     , {'H', 2.3, 'S',4.1  },
	# symmetry correction
	'sym'   , {'H',   0, 'S',-1.4 } );
my $infile = $ARGV[0] || die "# usage: $0 <promoters file>\n";
my %seqtable;
my $seq;
my $conttable = 0;
my $contseq = 0;
my $contador = 0;
open(SEQ, $infile) || die "# cannot open input $infile : $!\n";
while(<SEQ>)
{
	if(/^(b\d{4}) \\ ([ATGC]+)/)
	{
		#Se guardan los nombres y las secuencias en una estructura de datos (hash)
		my ($name,$seq) = ($1,$2); 
		$seqtable{$name} = $seq;
	}
}
close(SEQ);


# calculate NN free energy of a DNA duplex , dG(t) = (1000*dH - t*dS) / 1000
# parameters: 1) DNA sequence string; 2) Celsius temperature
# returns; 1) free energy scalar
# uses global hash %NNparams
my $total_total = 0;
my $totaldG = 0;
# En este ciclo se entra a la función duplex_deltaG en donde se realiza el conteo total de la energía libre de Gibbs que se encuentra en la secuencia en intervalos o ventanas de 15 nt. 
# También se calculan los valores de E1, E2 y D para predecir los promotores de acuerdo a sus niveles de corte.
foreach my $i(keys %seqtable){
		my $seq = $seqtable{$i};
		$totaldG = duplex_deltaG($seq, $T, $i);
		$total_total += $totaldG;
}

sub duplex_deltaG 
{
	#Variables utilizadas en la función 
	
	my ($seq,$tCelsius, $name) = @_;  #Parámetros pasados a la función desde el código principal.
	my ($DNAstep,$nt,$dG,$total_dG, $total_dG_partial) = ('','',0,0,0); #Inicialización de algunas variables
	my @sequence = split(//,uc($seq)); #Se corta la secuencia nucleótido por nucleótido y se guarda en un arreglo
	my @deltaG; # Arreglo en donde se guardarán los valores de la energía libre de Gibbs de las ventanas de 15 nt de las secuencias
	my $tK = 273.15 + $tCelsius;
	my $longitud = length($seq);
	my $k = 14;
	my $i = 0;
	my $mitad = 7;
	my $symetric = 1;
	sub complement{ $_[0] =~ tr/ATGC/TACG/; return $_[0] }	#FUnción que te regresa la cadena complementaria de DNA.
	# add dG for overlapping dinucleotides
	while ($k < length($seq)){
		$total_dG_partial = 0;
		for(my $n=$i;$n<$k;$n++) {
			$DNAstep = $sequence[$n].$sequence[$n+1].'/'.
				complement($sequence[$n].$sequence[$n+1]);
			if(!defined($NNparams{$DNAstep}))
			{
				$DNAstep = reverse($DNAstep);
			}
			$dG = ((1000*$NNparams{$DNAstep}{'H'})-
					($tK*$NNparams{$DNAstep}{'S'}))
					/ 1000 ;
			
			$total_dG_partial += $dG; 
		}
		# add correction for helix initiation
		$nt = $sequence[$i]; # first pair
		if(!defined($NNparams{$nt})){ $nt = complement($nt) } 
		$total_dG_partial += ((1000*$NNparams{$nt}{'H'})-
				($tK*$NNparams{$nt}{'S'}))
				/ 1000; 
		$nt = $sequence[$k]; # last pair
		if(!defined($NNparams{$nt})){ $nt = complement($nt) }
		$total_dG_partial += ((1000*$NNparams{$nt}{'H'})-
				($tK*$NNparams{$nt}{'S'}))
				/ 1000;	
		#Corrección por simetría
		my $a = $i;
		my $c = $k;
		while ($c > $mitad && $symetric == 1){
			my $complement = $sequence[$c];
			if($sequence[$a] eq complement($complement)){
				$symetric=1;
			}
			else {
				$symetric=0;
			}
			$c--;
			$a++;
		}
		if ($symetric == 1) {
			$total_dG_partial += ((1000*$NNparams{'sym'}{'H'})-
				($tK*$NNparams{'sym'}{'S'}))
				/ 1000; 	
		}
		# Se agregan los valores de la energía libre de Gibbs en un arreglo para utilizarlos posteriormente.
		push @deltaG, $total_dG_partial;
		$total_dG += $total_dG_partial;
		$symetric = 1;
		$i += 1;
		$k += 1;
		$mitad += 1;
	} 
	
	#Se calculan los valores de E1,E2 y D, para determinar, mediante valores de corte, si las secuencias pueden ser promotores. 
	my $min_1=0;
	my $max_1=50;
	my $min_2=100;
	my $max_2=200;
	my $x = 0;
	my $y= 0;
	my $contador1 = 0;
	while($max_2 <= scalar(@deltaG)){
		my $E1 = 0;
		my $E2 = 0;
		my $E1_partial = 0;
		my $E2_partial = 0;
		my $D = 0;
		my $min_real = 0;
		my $max_real = 0;
		#Suma de la energía libre de Gibbs de las primeras 50 ventanas
		for ($x=$min_1;$x<$max_1;$x++){
			$E1_partial+=$deltaG[$x];
		}
		#Suma de la energía libre de Gibbs de las 100 ventanas correspondientes
		for ($y=$min_2;$y<$max_2;$y++){
			$E2_partial+=$deltaG[$y];
		}
		$E1 = $E1_partial/50;
		$E2 = $E2_partial/100;
		$D = $E1 - $E2;
		if ($D > 2.76 && $E1 > -17.53){
			#Se ajustan los valores para que queden en términos de posiciones a nivel de secuencia de DNA
			$min_real = $max_1 - 400;
			$max_real = $min_2 - 400;
			print "$name \t $E1 \t $E2 \t $D \t $min_real \t $max_real\n";
		}
		$min_1+=1;
		$max_1+=1;
		$min_2+=1;
		$max_2+=1;
	}
	#Se regresa el valor total de la energía libre, solo para comprobar que los resultados estén bien.	
	return $total_dG;
}