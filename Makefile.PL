use 5.008000;
use ExtUtils::MakeMaker;

WriteMakefile(
	NAME		=> 'AI::NeuralNet::FastSOM',
	VERSION_FROM	=> 'lib/AI/NeuralNet/FastSOM.pm',
	ABSTRACT_FROM	=> 'lib/AI/NeuralNet/FastSOM.pm',
	AUTHOR		=> 'Rick Myers <jrm@cpan.org>',
	test => {
		TESTS => 't/*.t t/orig/*.t',
	},
	clean => {
		FILES => 't/*.bin',
	},
);

