use 5.006002;
use ExtUtils::MakeMaker;

WriteMakefile(
	NAME		=> 'AI::NeuralNet::FastSOM',
	VERSION_FROM	=> 'lib/AI/NeuralNet/FastSOM.pm',
	ABSTRACT_FROM	=> 'lib/AI/NeuralNet/FastSOM.pm',
	AUTHOR		=> 'Rick Myers <jrm@cpan.org>',
	PREREQ_PM	=> { Storable => 0 },
	TYPEMAPS	=> [ $] < 5.008000 ? 'typemap.v1' : 'typemap.v2' ],
	test		=> { TESTS => 't/*.t t/orig/*.t' },
	clean		=> { FILES => 't/*.bin' },
);

