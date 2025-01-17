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

my $n1 = $node_adaptor->fetch_by_taxon_id( 10 )
  || die "Could not retrieve node 10: $@";
ok( defined $n1,
	"Checking if the node 10 could be retrieved" );
is ($n1->taxon_id(),10,'Checking node ID is as expected');
diag($n1->to_string());

my $n2 = $node_adaptor->fetch_by_taxon_id( 8 )
  || die "Could not retrieve node 8: $@";
ok( defined $n2,
	"Checking if the node 8 could be retrieved" );
is ($n2->taxon_id(),8,'Checking node ID is as expected');
diag($n2->to_string());

ok($n1->has_ancestor($n2)==1,"Checking that ".$n1->to_string()." has ancestor ".$n2->to_string());
ok($n2->has_ancestor($n1)==0),"Checking that ".$n2->to_string()." does not have ancestor ".$n1->to_string();

done_testing;
