# Copyright [2009-2025] EMBL-European Bioinformatics Institute
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
									  
my $node_adaptor = Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyNodeAdaptor->new($tax);
		
ok( defined $node_adaptor, 'Checking if the node adaptor is defined' );

my $n1 = $node_adaptor->fetch_by_taxon_id( 10 )
  || die "Could not retrieve node 10: $@";
ok( defined $n1,
	"Checking if the node 10 could be retrieved" );
is ($n1->taxon_id(),10,'Checking node ID is as expected');
diag($n1->to_string());
$n1 = $node_adaptor->fetch_ancestor_by_rank($n1,"species");
diag($n1->to_string());

my $n2 = $node_adaptor->fetch_by_taxon_id( 11 )
  || die "Could not retrieve node 11: $@";
ok( defined $n2,
	"Checking if the node 11 could be retrieved" );
is ($n2->taxon_id(),11,'Checking node ID is as expected');
diag($n2->to_string());
$n2 = $node_adaptor->fetch_ancestor_by_rank($n2,"species");
diag($n2->to_string());

done_testing;
