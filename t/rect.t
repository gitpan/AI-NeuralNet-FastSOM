#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw(no_plan);
BEGIN { use_ok('AI::NeuralNet::FastSOM::Rect') };

######
#use Data::Dumper;

{
	my $nn = AI::NeuralNet::FastSOM::Rect->new(
		output_dim => '5x6',
		input_dim  => 3
	);

	ok( $nn->isa( 'AI::NeuralNet::FastSOM::Rect' ), 'rect class' );

	my $nn2 = $nn;
	my $nn3 = $nn2;
	is( $nn, $nn3, 'rect eq' );
	is( $nn->refcount(), 3, 'rect refs' );

	my $m1 = $nn->map;
	isa_ok( $m1, 'ARRAY', 'map array' );

	my $m2 = $m1;
	my $m3 = $nn2->map;
	my $m4 = $m3;
	is( $m2, $m4, 'map eq' );

	my $a = $m1->[0];
	isa_ok( $a, 'ARRAY', 'array array' );
	ok( $a != $m1, 'array unique' );

	my $a2 = $m4->[0];
	if ( !is( $a, $a2, 'array eq' ) ) {
		warn "not equal:\n";
		warn "    $a\n    $a2\n";
	}

	my $v = $a->[0];
	isa_ok( $v, 'ARRAY', 'vector array' );
	ok( $v != $a, 'vector unique' );

	my $v2 = $nn3->map->[0]->[0];
	if ( !is( $v, $v2, 'vector eq' ) ) {
		warn "not equal:\n";
		warn "    $v\n";
		warn "    $v2\n";
	}
	my $v3 = $nn2->map->[0][0];
	if ( !is( $v, $v3, 'vector shorter' ) ) {
		warn "not equal:\n";
		warn "    $v\n    $v3\n";
	}

	my $m = $nn->map;
	$m->[0][0][0] = 3.245;
	is( $m->[0][0][0], 3.245, 'element set' );
	$m->[0][0][0] = 1.25;
	is( $m->[0][0][0], 1.25, 'element reset' );
	$m->[0][0][1] = 4.8;
	is( $m->[0][0][1], 4.8, 'element set z' );
	$m->[0][0][1] = 2.6;
	is( $m->[0][0][1], 2.6, 'element reset z' );
	$m->[0][1][0] = 8.9;
	is( $m->[0][1][0], 8.9, 'element set y' );
	$m->[0][1][0] = 1.2;
	is( $m->[0][1][0], 1.2, 'element reset y' );
	$m->[1][0][0] = 5.4;
	is( $m->[1][0][0], 5.4, 'element set z' );
	$m->[1][0][0] = 3.23;
	is( $m->[1][0][0], 3.23, 'element reset z');

	$m->[4][5][2] = 2.29;
	is( $m->[4][5][2], 2.29, 'last element set' );
	is( $m->[-1][5][2], 2.29, 'negative x' );
	is( $m->[4][-1][2], 2.29, 'negative y' );
	is( $m->[4][5][-1], 2.29, 'negative z' );
	is( $m->[-1][-1][-1], 2.29, 'negative all' );

}

{
	my $nn = AI::NeuralNet::FastSOM::Rect->new(
		output_dim => '5x6',
		input_dim  => 3
	);
	ok ($nn->isa ('AI::NeuralNet::FastSOM::Rect'), 'class');
	is ($nn->{_X}, 5, 'X');
	is ($nn->{_Y}, 6, 'Y');
	is ($nn->{_Z}, 3, 'Z');
	is ($nn->radius, 2.5, 'radius');
	is ($nn->output_dim, "5x6", 'output dim');
}

{
	my $nn = new AI::NeuralNet::FastSOM::Rect(
		output_dim => "5x6",
		input_dim  => 3
	);
	$nn->initialize;
#	print Dumper $nn;
#	exit;

	my @vs = ([ 3, 2, 4 ], [ -1, -1, -1 ], [ 0, 4, -3]);
	$nn->train(400, @vs);

	foreach my $v (@vs) {
		ok(_find($v,$nn->map),'found learned vector '.join(",", @$v));
	}

	sub _find {
		my $v = shift;
		my $m = shift;

		use AI::NeuralNet::FastSOM::Utils;
		foreach my $x ( 0 .. 4 ) {
			foreach my $y ( 0 .. 5 ) {
				my $rv = AI::NeuralNet::FastSOM::Utils::vector_distance($m->[$x]->[$y], $v);
				return 1 if $rv < 0.01;
			}
		}
		return 0;
	}

	ok ($nn->as_string, 'pretty print');
	ok ($nn->as_data, 'raw format');

#	print $nn->as_string;
}

{
    my $nn = new AI::NeuralNet::FastSOM::Rect (output_dim => "5x6",
					   input_dim  => 3);
    $nn->initialize;

    foreach my $x (0 .. 5 -1) {
	foreach my $y (0 .. 6 -1 ) {
	    ok ( (!grep { $_ > 0.5 || $_ < -0.5 } @{ $nn->value ( $x, $y ) }) , "$x, $y: random vectors in [-0.5, 0.5]");
	}
    }
}

{
	my $nn = new AI::NeuralNet::FastSOM::Rect(
		output_dim => "5x6",
		input_dim  => 3
	);
	$nn->initialize;
#       print Dumper $nn;
#       exit;

	my @vs = ([ 3, 2, 4 ], [ -1, -1, -1 ], [ 0, 4, -3]);
	$nn->train(400, @vs);

	my $k = keys %$nn;
	is( $k, 10, 'scalar rect key count' );
	my @k = keys %$nn;
	is( @k, 10, 'array rect key count' );
}

use File::Temp 'tempfile';
use Storable qw/store retrieve/;
my ($file,$bmu_x,$bmu_y);
{
	my $nn = AI::NeuralNet::FastSOM::Rect->new(
		output_dim => '5x6',
		input_dim  => 3
	);
	$nn->initialize;

	my @vs = ([ 3, 2, 4 ], [ -1, -1, -1 ], [ 0, 4, -3]);
	$nn->train(400, @vs);

	($bmu_x,$bmu_y) = $nn->bmu([3,2,4]);

	$file = tempfile();
	store( $nn, $file );
}
{
	my $nn = retrieve( $file );
	my ($x,$y) = $nn->bmu([3,2,4]);
	is( $x, $bmu_x, 'stored x' );
	is( $y, $bmu_y, 'stored y' );
	unlink $file;
}

__END__

# randomized pick
    @vectors = ...;
my $get = sub {
    return @vectors [ int (rand (scalar @vectors) ) ];
    
}
$nn->train ($get);

# take exactly 500, round robin, in order
our $i = 0;
my $get = sub {
    return undef unless $i < 500;
return @vectors [ $i++ % scalar @vectors ];
}
