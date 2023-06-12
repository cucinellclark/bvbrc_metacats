#The MSA application with variance analysis.

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;

use strict;
use P3DataAPI;
use Data::Dumper;
use File::Basename;
use File::Slurp;
use LWP::UserAgent;
use JSON::XS;
use JSON;
use IPC::Run qw(run);
use Cwd;
use Clone;
use URI::Escape;

my $script = Bio::KBase::AppService::AppScript->new(\&process_metacats);
my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;

my $rc = $script->run(\@ARGV);

exit $rc;


sub process_metacats
{
    my($app, $app_def, $raw_params, $params) = @_;

    print "Proc Metacats ", Dumper($app_def, $raw_params, $params);
    my $token = $app->token();
    my $data_api_module = P3DataAPI->new($data_api, $token);
    my $output_folder = $app->result_folder();

    #
    # Create an output directory under the current dir. App service is meant to invoke
    # the app script in a working directory; we create a folder here to encapsulate
    # the job output.
    #
    # We also create a staging directory for the input files from the workspace.
    #

    my $cwd = getcwd();
    my $work_dir = "$cwd/work";
    my $stage_dir = "$cwd/stage";

    -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
    -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";

    my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
    my $dat = { data_api => $data_api };
    my $sstring = encode_json($dat);

    #
    # Read parameters and discover input files that need to be staged.
    #
    # Make a clone so we can maintain a list of refs to the paths to be
    # rewritten.
    #
    my %in_files;
    my $params_to_app = Clone::clone($params);
    # Get values common to all inputs.
    my $prefix = $params_to_app->{output_file};
    my $p_value = $params_to_app->{p_value};
    my $alphabet = $params_to_app->{alphabet};
    my $input_type = $params_to_app->{input_type};
    my $check_header = 0;
    my $seqFile = "";
    my $metaDataFile = "";
    my $staged = {};
    my $ofile = "$stage_dir/downloaded_sequences.fasta";

    # Get values based on input type.
    if ($input_type eq "files") {
        # Get sequence file and metadata file in the staging directory.
        my @to_stage;
        $check_header = 1;
        push(@to_stage, $params_to_app->{alignment_file});
        push(@to_stage, $params_to_app->{group_file});
        if (@to_stage)
        {
            warn Dumper(\@to_stage);
            $staged = $app->stage_in(\@to_stage, $stage_dir, 1);
            my %new_hash = %{ $staged };
            $seqFile = $new_hash{$params_to_app->{alignment_file}};
            $metaDataFile = $new_hash{$params_to_app->{group_file}};
        }
    } elsif ($input_type eq "groups") {
        # Get the sequences and create a metadata file for the groups.
        open(F, ">$ofile") or die "Could not open $ofile";
        $metaDataFile = "$work_dir/metadata.tsv";
        open(G, ">$metaDataFile") or die "Could not open $metaDataFile";
        print G "Seq_ID\t$prefix\n";
        for my $feature_name (@{$params_to_app->{groups}}) {
            print STDOUT "Getting features in $feature_name\n";
            my $ids = $data_api_module->retrieve_patricids_from_feature_group($feature_name);
            my $seq = "";
            if ($alphabet eq "na") {
                $seq = $data_api_module->retrieve_nucleotide_feature_sequence($ids);
            } else {
                $seq = $data_api_module->retrieve_protein_feature_sequence($ids);
            }
            for my $id (@$ids) {
                my $out = ">$id\n" . uc($seq->{$id}) . "\n";
                print F $out;
                print G "$id\t" . "$feature_name" . "\n";
            }
        }
        close F;
        close G;
    } elsif ($input_type eq "auto") {
        # Put the sequences in the ofile from the patric ids, and get the metadata file from the JSON object.
        open(F, ">$ofile") or die "Could not open $ofile";
        $metaDataFile = "$work_dir/metadata.tsv";
        open(G, ">$metaDataFile") or die "Could not open $metaDataFile";
        print G "Seq_ID\t$prefix\n";
        my @ids = ();
        for (@{$params->{auto_groups}}) {
            push(@ids, $_->{id});
            print G $_->{id} . "\t" . $_->{grp} . "\n";
        }
        my $seq = "";
        if ($alphabet eq "na") {
            $seq = $data_api_module->retrieve_nucleotide_feature_sequence(\@ids);
        } else {
            $seq = $data_api_module->retrieve_protein_feature_sequence(\@ids);
        }
        for my $id (@ids) {
            print F ">$id\n" . uc($seq->{$id}) . "\n";
        }
        close F;
        close G;
    } elsif ($input_type eq "list") {
        # Put the sequences in the ofile from the patric ids, and get the metadata file from the JSON object.
        open (F, ">$ofile") or die "Couldnot open $ofile";
        $metaDataFile = "$work_dir/metadata.tsv";
        open(G, ">$metaDataFile") or die "Could not open $metaDataFile";
        print G "Seq_ID\t$prefix\n";
        my $seq = "";
        my @ids = @{$params->{feature_list}};
        if ($alphabet eq "na") {
            $seq = $data_api_module->retrieve_nucleotide_feature_sequence(\@ids);
        } else {
            $seq = $data_api_module->retrieve_protein_feature_sequence(\@ids);
        }
        for my $id (@ids) {
            print F ">$id\n" . uc($seq->{$id}) . "\n";
        }
        close F;
        close G;
    } else {
        die("Unrecognized input type.");
    }
    if (($input_type eq "groups") or ($input_type eq "auto") or ($input_type eq "list")) {
        # Align the sequences.
        my @mafft_cmd = ("mafft", "--auto", "--preservecase", $ofile);
        my $string_cmd = join(" ", @mafft_cmd);
        print STDOUT "Running mafft.\n";
        print STDOUT "$string_cmd\n";
        $seqFile = "$work_dir/output.afa";
        my $ok = run(\@mafft_cmd, "1>", $seqFile, "2>", "$work_dir/$prefix.mafft.log");
        if (!$ok) {
            die "Mafft command failed.\n";
        }
        print STDOUT "Finished mafft.\n";
    }
    # Replace leading and trailing gaps with the '#' symbol in the sequence file.
    open my $fh, '<', "$seqFile" or die "Cannot open $seqFile: $!";
    my $adjusted_seqFile = "$work_dir/adjusted_seqs.fasta";
    open(OUT, '>', "$adjusted_seqFile") or die "Cannot open $adjusted_seqFile: $!";
    my $seq_string = "";
    my $header = "";
    while ( my $line = <$fh> ) {
        chomp $line;
        if (substr($line, 0, 1) eq ">") {
            if ($header) {
                print OUT "$header\n";
            }
            $header = $line;
            if ($seq_string) {
                $seq_string = replace_gaps($seq_string);
                print OUT "$seq_string\n";
            }
            $seq_string = "";
        } else {
            $seq_string = $seq_string.$line;
        }
    }
    if ($header) {
        print OUT "$header\n";
    }
    if ($seq_string) {
        $seq_string = replace_gaps($seq_string);
        print OUT "$seq_string\n";
    }
    close OUT;
    # Run the analysis.
    my @cmd = ("metadata_parser", $adjusted_seqFile, $metaDataFile, $alphabet, $p_value, $check_header,  "$work_dir/");
    run_cmd(\@cmd);
    my @output_suffixes = (
        [qr/Table\.tsv$/, "tsv"],
        [qr/\.log$/, "txt"],
        [qr/\.afa$/, "aligned_protein_fasta"]
        );
    opendir(D, $work_dir) or die "Cannot opendir $work_dir: $!";
    my @files = sort { $a cmp $b } grep { -f "$work_dir/$_" } readdir(D);
    my $output=1;
    for my $file (@files)
    {
	for my $suf (@output_suffixes)
	{
	    if ($file =~ $suf->[0])
	    {
 	    	$output=0;
		my $path = "$output_folder/$file";
		my $type = $suf->[1];
		$app->workspace->save_file_to_file("$work_dir/$file", {}, "$output_folder/$file", $type, 1,
					       (-s "$work_dir/$file" > 10_000 ? 1 : 0), # use shock for larger files
					       $token);
	    }
	}
    }
    #
    # Clean up staged input files.
    #
    while (my($orig, $staged_file) = each %$staged)
    {
	unlink($staged_file) or warn "Unable to unlink $staged_file: $!";
    }
    return $output;
}

sub replace_gaps {
    my($line) = @_;
    my $start_len = 0;
    my $end_len = 0;
    if ($line =~ /^(-+)/) {
        $start_len = length($1);
    }
    if ($line =~ /(-+)$/) {
       $end_len = length($1);
    }
    substr($line, 0, $start_len) = '#' x $start_len;
    substr($line, length($line) - $end_len, $end_len) = '#' x $end_len;
    return $line;
}

sub run_cmd() {
    my $cmd = $_[0];
    my $ok = run(@$cmd);
    if (!$ok)
    {
        die "Command failed: @$cmd\n";
    }
}
