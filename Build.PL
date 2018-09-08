use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'MTToolbox',
    license             => 'perl',
    dist_author         => 'Scott Yourstone <scott.yourstone81@gmail.com>',
    dist_version_from   => 'lib/MTToolbox.pm',
    requires => {
        'Test::More' => 0,
        'version'    => 0,
    },
    build_requires => {
		'Carp' => 0,
		'Class::Std::Utils' => 0,
		'Chart::Gnuplot' => 0,
		'Chart::Graph::Gnuplot' => 0,
		'Data::Dumper' => 0,
		'Exception::Class' => 0,
		'File::Basename' => 0,
		'File::Spec' => 0,
		'File::Temp' => 0,
		'Getopt::Long' => 0,
		'IPC::Cmd' => 0,
		'List::MoreUtils' => 0,
		'Math::BigInt' => 0,
		'Module::Build' => 0,
		'Pod::Coverage' => 0,
		'Scalar::Util' => 0,
		'Statistics::Descriptive' => 0,
		'Test::Exception' => 0,
		'Test::More' => 0,
		'Test::Pod' => 0,
		'Test::Pod::Coverage' => 0,
		'Test::Warn' => 0,
		'Tk' => 0,
		'Tk::DialogBox' => 0,
		'Tk::LabEntry' => 0,
		'Tk::NoteBook' => 0,
		'XML::Parser' => 0,
		'XML::Simple' => 0,
		'version' => 0,
    },
    add_to_cleanup      => [ 'MTToolbox-*' ],
    recursive_test_files => 1,
);

$builder->create_build_script();
