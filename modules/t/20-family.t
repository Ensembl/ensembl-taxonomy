# Copyright [2009-2023] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Taxonomy::TaxonomyNode;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyNodeAdaptor;
use Bio::EnsEMBL::Test::MultiTestDB;

use FindBin qw($Bin);

my $conf_file = "$Bin/db.conf";
my $conf = do $conf_file
  || die "Could not load configuration from " . $conf_file;

$conf = $conf->{tax_test};

my $test_db_dir = $FindBin::Bin;
my $testdb  = Bio::EnsEMBL::Test::MultiTestDB->new('multi', $test_db_dir);

my $dba = $testdb->get_DBAdaptor('taxonomy');

my $node_adaptor = Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyNodeAdaptor->new($dba);

ok( defined $node_adaptor, 'Checking if the node adaptor is defined' );

my $node = $node_adaptor->fetch_by_taxon_id( $conf->{taxon_id} )
  || die "Could not retrieve node $conf->{taxon_id}: $@";
ok( defined $node,
	"Checking if the node node $conf->{taxon_id} could be retrieved" );
is ($node->taxon_id(),$conf->{taxon_id},'Checking node ID is as expected');

# get lineage
my @lineage = @{$node_adaptor->fetch_ancestors($node)};
ok ((scalar @lineage)>1,'Checking lineage is found');

# get descendants
my @desc = @{$node_adaptor->fetch_descendants($node)};
ok ((scalar @desc)>0,'Checking descendants are found');
# get leaf nodes
my @leaves = @{$node_adaptor->fetch_leaf_nodes($node)};
ok ((scalar @leaves)>0,'Checking leaves are found');

done_testing;
