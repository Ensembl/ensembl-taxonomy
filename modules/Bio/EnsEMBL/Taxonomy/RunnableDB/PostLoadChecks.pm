=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Taxonomy::RunnableDB::PostLoadChecks

=head1 DESCRIPTION

Sanity checks on the final NCBI taxonomy database. 
First check to make sure the total count of taxa nodes match the sum of the left and right indexes
Second check find duplicated items in the ncbi_taxa_name table, report and clean them up

=cut

package Bio::EnsEMBL::Taxonomy::RunnableDB::PostLoadChecks;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Production::Pipeline::Common::Base/;

sub fetch_input {
    my $self = shift;

}

sub run {
    my $self = shift;
    # Checking that the total count is matching the sum of left or right indexes
    my $sth_total_count = $self->count_rows("SELECT COUNT(*) FROM ncbi_taxa_node");
    my $sth_left_index = $self->count_rows("SELECT SUM(left_index > 0) FROM ncbi_taxa_node");
    my $sth_right_index = $self->count_rows("SELECT SUM(right_index > 0) FROM ncbi_taxa_node");
    # If all the counts are the same, it's fine.
    if ($sth_total_count==$sth_left_index and $sth_total_count==$sth_right_index and $sth_left_index==$sth_right_index)
    {
      $self->warning('Count '.$sth_total_count.' is matching left indexes sum '.$sth_left_index.' and right indexes sum '.$sth_right_index);
    }
    else
    {
       $self->throw('Count '.$sth_total_count.' is not matching left indexes sum '.$sth_left_index.' and right indexes sum '.$sth_right_index.' .Please contact the Compara team');
    }
    #Checking that there is no duplicated rows in the ncbi_taxa_name table of the ncbi taxonomy database.
    my $sql_duplicates = 'SELECT taxon_id, name, name_class, count(*) AS duplicate_count FROM ncbi_taxa_name GROUP BY taxon_id, name, name_class HAVING duplicate_count > 1;';
    my $sth_duplicates = $self->data_dbc->prepare($sql_duplicates);
    $sth_duplicates->execute();
    #If there are any duplicated rows, write a message in the hive msg table and delete the duplicates from the ncbi_taxa_name table.
    #if ($sth_duplicates->fetchrow_array())
    #{
    my $duplicates = 0;

    # Getting all the duplicate items
    my $all_duplicates = $sth_duplicates->fetchall_arrayref({});
    foreach my $duplicate (@$all_duplicates) {
      # Calculate number of rows to delete
      my $number_of_row_to_delete=$duplicate->{"duplicate_count"} - 1;
      # Write in the hive msg table
      $self->warning('The taxon id '.$duplicate->{"taxon_id"}.' had '.$number_of_row_to_delete.' duplicates for "'.$duplicate->{"name_class"}.'" named "'.$duplicate->{"name"}.'". The duplicates have been cleaned up in the database.');
      # Delete duplicates
      my $sql_delete_duplicates= 'DELETE from ncbi_taxa_name where name=? and taxon_id=? and name_class=?  LIMIT ?;';
      my $sth_delete_duplicates = $self->data_dbc->prepare($sql_delete_duplicates);
      $sth_delete_duplicates->execute($duplicate->{"name"},$duplicate->{"taxon_id"},$duplicate->{"name_class"},$number_of_row_to_delete);

      $duplicates = 1;
    }

    if ($duplicates) {
      # Look for duplicates again, to check that all the fixes worked.
      $sth_duplicates->execute();
      if ($sth_duplicates->fetchrow_array()) {
        $self->throw('Duplicate names remain after trying to remove them.');
      }
    } else {
      # Write in the hive msg table if there is no duplicate
      $self->warning('Found no duplicates in the database')
    }
}

# subroutine to run a count query and return the result
sub count_rows {

  my ($self, $sql) = @_;
  my $sth = $self->data_dbc->prepare($sql);
  $sth->execute();

  return ($sth->fetchrow_array());

}

1;
