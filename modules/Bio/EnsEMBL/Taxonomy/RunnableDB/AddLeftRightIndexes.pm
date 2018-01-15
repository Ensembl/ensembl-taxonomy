=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=pod 

=head1 NAME

Bio::EnsEMBL::Taxonomy::RunnableDB::AddLeftRightIndexes

=head1 DESCRIPTION

This Runnable adds left_index and right_index to an NCBI database.
It modifies the database listed in the 'db_conn' parameter or the current hive database.

=cut

package Bio::EnsEMBL::Taxonomy::RunnableDB::AddLeftRightIndexes;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');

# Runnable interface
#####################

sub param_defaults {
    return {
        'root_node_name'    => 'root',      # The name of the root node in the taxonomy. We'll start indexing from there. It could in fact be any node
    }
}

sub fetch_input {
    my $self = shift;

    # Stetement to get the children of a given node
    my $sql_get_children = 'SELECT taxon_id FROM ncbi_taxa_node WHERE parent_id = ?';
    my $sth_get_children = $self->data_dbc->prepare($sql_get_children);
    $self->param('sth_get_children', $sth_get_children);

    # Statement to update left/right_index of a node
    my $sql_update_node = 'UPDATE ncbi_taxa_node SET left_index = ?, right_index = ? WHERE taxon_id = ?';
    my $sth_update_node = $self->data_dbc->prepare($sql_update_node);
    $self->param('sth_update_node', $sth_update_node);

    # Get the root node
    my $sql_get_root_node = 'SELECT ncbi_taxa_node.* FROM ncbi_taxa_node JOIN ncbi_taxa_name USING (taxon_id) WHERE ncbi_taxa_name.name = ? LIMIT 1';
    my $root_node = $self->data_dbc->db_handle->selectrow_hashref($sql_get_root_node, undef, $self->param_required('root_node_name'));
    die sprintf("Unable to find a node with name '%s' in '%s'\n", $self->param('root_node_name'), $self->data_dbc->locator) unless $root_node;
    $self->param('root_node', $root_node);
}

sub run {
    my $self = shift;
    $self->build_store_leftright_indexing($self->param('root_node'), 1);  # The first left_index is 1
}


# Internal methods
###################

sub build_store_leftright_indexing {
    my ($self, $node, $counter) = @_;

    # Get all the children
    # NOTE: since the same statement handle is going to be used throughout
    # the recursion, we need to pull all the data and "finish" it before
    # recursing, hence the "fetchall" here instead of a while-fetch-recurs
    # loop below
    my $sth_get_children = $self->param('sth_get_children');
    $sth_get_children->execute($node->{'taxon_id'});
    my $all_nodes = $sth_get_children->fetchall_arrayref({});
    $sth_get_children->finish();

    my $left_index = $counter++;
    foreach my $child_node (@$all_nodes) {
        $counter = $self->build_store_leftright_indexing($child_node, $counter);
    }
    my $right_index = $counter++;

    # Store the new values;
    warn sprintf('taxon_id %d: %d-%d', $node->{'taxon_id'}, $left_index, $right_index), "\n" if $self->debug;
    $self->param('sth_update_node')->execute($left_index, $right_index, $node->{'taxon_id'});

    return $counter;
}

1;
