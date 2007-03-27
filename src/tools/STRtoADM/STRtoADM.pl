## Program to convert from the Structure format to the Admixmap format for the genotype data.

## By Indrani Halder (iindrani@hotmail.com)

## run with perl STRtoADM.pl <inputfile> [<outputfile>]
## default outputfile name is "genotypes.txt"
## Intermediate file called output.txt is also created.

## If your genotypes file is not coded in numbers then use a text editor to convert alleles to numbers first.
## Code missing data as "?" , for both alleles at a locus. See example inputfile.

$inputfile =  $ARGV[0];#"inputfile.txt";
$outputfile = "genotypes.txt";
if( $#ARGV > 0){$outputfile =  $ARGV[1];}

open(FILE, $inputfile);
open(OFILE,">output.txt");
$data=<FILE>;
chomp($data);
print OFILE $data,"\n";
while($data=<FILE>)
{
	chomp($data);
	@da=split(/\s+/,$data);
	$len=@da;
	print OFILE $da[0],"\t";
	for($i=1;$i<$len;$i+=2)
	{
		print OFILE $da[$i],$da[$i+1],"\t";
	}
	print OFILE "\n";
}
close(OFILE);
close(FILE);

# this generates the intermediate file
# now the intermediate file is converted to Admixmap file

open(FILE,"output.txt");

open(OFILE,">$outputfile");

$data=<FILE>;
chomp($data);
@da = split(/\t/,$data);
$len = @da;

for($i=0;$i<$len-1;$i++)
{
	print OFILE "\"",$da[$i],"\" ";## print locus names, in quotes, followed by space
}
print OFILE "\"",$da[$len-1],"\"";## print final locus name without space
print OFILE "\n";

while($data = <FILE>)
{
	chomp($data);
	@da = split(/\t/,$data);
	$len = @da;
	print OFILE "\"",$da[0],"\" ";## print Individual id, in quotes
	for($i=1;$i<$len;$i++)
	{
		$v1=substr($da[$i],0,1);## allele1
		$v2=substr($da[$i],1,1);## allele2
		if($v1 ne "?")
		{
			print OFILE "\"",$v1,",",$v2,"\"";
		}
		else## missing genotype
		{
			print OFILE "\"\"";
		}
		if($i<$len-1)##print space if not last locus
		  {
		    print OFILE " ";
		  }

	}
	print OFILE "\n";
}

close(OFILE);
close(FILE);
