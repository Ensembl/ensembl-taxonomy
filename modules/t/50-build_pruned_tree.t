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

my @expected;
my @unexpected;
my $sub = sub {
	my ($node,$depth) = @_;
	my $found=0;
	for(my $i=0; $i<scalar @expected; $i++) {
		if(defined $expected[$i] && $expected[$i]==$node->taxon_id()) {
			delete $expected[$i];
			$found=1;
			last;
		}
	}
	if($found==0) {
		push @unexpected,$node->taxon_id();
	}
	return;
};

my $root = $node_adaptor->fetch_by_taxon_id(1);
$node_adaptor->add_descendants($root);
diag("Total tree:", $root->to_tree_string());

# test 1		
my $leaves = $node_adaptor->fetch_all_by_taxon_ids([10]);
$root = $node_adaptor->fetch_by_taxon_id(1);
diag("Tree before [1=>10]:", $root->to_tree_string());
ok($root->taxon_id()==1,"Checking root id");
@expected = (1,2,4,8,10);
diag("Expected: ",join(' ',grep{defined}@expected));
$node_adaptor->build_pruned_tree($root,$leaves);
diag("Tree after:", $root->to_tree_string());
@unexpected = ();
$root->traverse_tree($sub);
is(scalar @expected,0,"All expected found");
diag("Expected not found: ",join(' ',grep{defined}@expected));
is(scalar @unexpected,0,"No unexpected found");
diag("Unexpected found: ",join(' ',@unexpected));

# test 2
$leaves = $node_adaptor->fetch_all_by_taxon_ids([10,9,6]);
$root = $node_adaptor->fetch_by_taxon_id(1);
ok($root->taxon_id()==1,"Checking root id");
diag("Tree before [1=>10,9,6]:", $root->to_tree_string());
@expected = (1,2,4,8,10,3,5,9,6);
diag("Expected: ",join(' ',grep{defined}@expected));
$node_adaptor->build_pruned_tree($root,$leaves);
diag("Tree after:\n", $root->to_tree_string());
@unexpected = ();
$root->traverse_tree($sub);
is(scalar @expected,0,"All expected found");
diag("Expected not found: ",join(' ',grep{defined}@expected));
is(scalar @unexpected,0,"No unexpected found");
diag("Unexpected found: ",join(' ',@unexpected));

# test 3
$leaves = $node_adaptor->fetch_all_by_taxon_ids([7,8]);
$root = $node_adaptor->fetch_by_taxon_id(2);
ok($root->taxon_id()==2,"Checking root id");
diag("Tree before [2=>7,8]:", $root->to_tree_string());
@expected = (2,4,8,7);
diag("Expected: ",join(' ',grep{defined}@expected));
$node_adaptor->build_pruned_tree($root,$leaves);
diag("Tree after:\n", $root->to_tree_string());
@unexpected = ();
$root->traverse_tree($sub);
is(scalar @expected,0,"All expected found");
diag("Expected not found: ",join(' ',grep{defined}@expected));
is(scalar @unexpected,0,"No unexpected found");
diag("Unexpected found: ",join(' ',@unexpected));

# test 4		
$leaves = $node_adaptor->fetch_all_by_taxon_ids([8,10]);
$root = $node_adaptor->fetch_by_taxon_id(1);
ok($root->taxon_id()==1,"Checking root id");
diag("Tree before [1=>8,10]:", $root->to_tree_string());
@expected = (1,2,4,8,10);
diag("Expected: ",join(' ',grep{defined}@expected));
$node_adaptor->build_pruned_tree($root,$leaves);
diag("Tree after:\n", $root->to_tree_string());
@unexpected = ();
$root->traverse_tree($sub);
is(scalar @expected,0,"All expected found");
diag("Expected not found: ",join(' ',grep{defined}@expected));
is(scalar @unexpected,0,"No unexpected found");
diag("Unexpected found: ",join(' ',@unexpected));

# test 5		
$leaves = $node_adaptor->fetch_all_by_taxon_ids([8,10,9,11]);
for my $leaf (@{$leaves}) {
	$leaf->{dba} = [1];
}
$root = $node_adaptor->fetch_by_taxon_id(1);
ok($root->taxon_id()==1,"Checking root id");
diag("Tree before [1=>8,9,10,11]:", $root->to_tree_string());
@expected = (1,2,3,4,8,10,5,9,11);
diag("Expected: ",join(' ',grep{defined}@expected));
$node_adaptor->build_pruned_tree($root,$leaves);
diag("Tree after:\n", $root->to_tree_string());
@unexpected = ();
$root->traverse_tree($sub);
is(scalar @expected,0,"All expected found");
diag("Expected not found: ",join(' ',grep{defined}@expected));
is(scalar @unexpected,0,"No unexpected found");
diag("Unexpected found: ",join(' ',@unexpected));
$node_adaptor->collapse_tree($root);
diag("Tree after collapse:\n", $root->to_tree_string());

done_testing;


