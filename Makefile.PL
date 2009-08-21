use 5.006002;
use ExtUtils::MakeMaker;

WriteMakefile(
	NAME		=> 'AI::NeuralNet::FastSOM',
	VERSION_FROM	=> 'lib/AI/NeuralNet/FastSOM.pm',
	ABSTRACT_FROM	=> 'lib/AI/NeuralNet/FastSOM.pm',
	AUTHOR		=> 'Rick Myers <jrm@cpan.org>',
	PREREQ_PM	=> { Storable => 0 },
	#TYPEMAPS	=> [ $] < 5.008000 ? 'typemap.v1' : 'typemap.v2' ],
	test		=> { TESTS => 't/*.t t/orig/*.t' },
	clean		=> { FILES => 't/*.bin typemap' },
);

#
# everything below is a work-around for some sort of bug in ExtUtils::ParseXS
# not picking up typemap files unless named "typemap" in perl5.6.2
#

package MY;

sub xs_c {
	my $t = shift->SUPER::xs_c(@_);
	$t =~ s/:/:\n	\$(MAKE) typemap/;
	$t;
}

sub postamble {
	if ( $] < 5.008000 ) {
	return <<'EOP';
typemap:
	$(CP) typemap.v1 typemap
EOP
	}
	return <<'EOP';
typemap:
	$(CP) typemap.v2 typemap
EOP
}

