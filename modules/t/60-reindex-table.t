# Copyright [2009-2014] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Bio::EnsEMBL::Hive::Utils::Test qw(standaloneJob);

my $root_node_name = 'Bacteria';

plan tests => 4;

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('eg');
my $dba = $multi_db->get_DBAdaptor('tax');
my $dbc = Bio::EnsEMBL::Hive::DBSQL::DBConnection->new(-dbconn => $dba->dbc);

# We reset the indices
$dbc->do('UPDATE ncbi_taxa_node SET left_index = 0, right_index = 0');

# And we rebuild them
standaloneJob(
    'Bio::EnsEMBL::RunnableDB::AddLeftRightIndexes',
    {
        'db_conn'   => $dbc->url,       # Parameters in the test still have to be stringified to mimic the job.input_id table
        'root_node_name'    => $root_node_name, # Searches "root" otherwise
    },
);

my $node_adaptor = Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor->new($dba);
ok( defined $node_adaptor, 'Checking if the node adaptor is defined' );

my $root = $node_adaptor->fetch_all_by_name($root_node_name)->[0];
ok( defined $root, 'Could still find the root in the database' );

my $sub = sub {
    my ($node, $depth) = @_;

    ok($node->{'left_index'} < $node->{'right_index'}, 'left_index is smaller than right_index');

    my $index_diff = $node->{'right_index'} - $node->{'left_index'};
    ok($index_diff % 2, 'The difference between left_index and right_index is odd');

    my $desc_ids = $node_adaptor->fetch_descendant_ids($node);
    my $expected_num_descendants = ($index_diff - 1) / 2;   # guaranteed to be an integer
    is($expected_num_descendants, scalar(@$desc_ids), 'There is the correct number of descendants');
};

my $all_descendants = $node_adaptor->fetch_descendant_ids($root);
subtest 'left_index and right_index are correctly set' => sub
{
    plan tests => 3*(scalar(@$all_descendants)+1);
    $root->traverse_tree($sub);
};

done_testing();

