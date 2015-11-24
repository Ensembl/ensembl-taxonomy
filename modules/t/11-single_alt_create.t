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

use Bio::EnsEMBL::Taxonomy::TaxonomyNode;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyNodeAdaptor;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('multi');
my $tax   = $multi->get_DBAdaptor('taxonomy');

use FindBin qw($Bin);
my $conf_file = "$Bin/db.conf";

my $conf = do $conf_file
  || die "Could not load configuration from " . $conf_file;
  
$conf = $conf->{tax_test};
									  
my $node_adaptor = $tax->get_TaxonomyNodeAdaptor();
		
ok( defined $node_adaptor, 'Checking if the node adaptor is defined' );

my $node = $node_adaptor->fetch_by_taxon_id( $conf->{taxon_id} )
  || die "Could not retrieve node $conf->{taxon_id}: $@";
ok( defined $node,
	"Checking if the node node $conf->{taxon_id} could be retrieved" );
is ($node->taxon_id(),$conf->{taxon_id},'Checking node ID is as expected');
ok ((scalar keys %{$node->names()})>0,'Checking node names are populated');
ok (defined $node->rank(),'Checking node rank is set');

# now try and find parents and children
my $parent = $node->parent();
ok (defined $parent,'Checking parent is found');

my $children = $node->children();
ok ((scalar @$children)>0,'Checking children are found');
done_testing;
