=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 NAME

Bio::EnsEMBL::Taxonomy::PipeConfig::ImportNCBItaxonomy_conf

=head1 DESCRIPTION

A pipeline to import NCBI taxonomy data into a ncbi_taxonomy database

=cut

package Bio::EnsEMBL::Taxonomy::PipeConfig::ImportNCBItaxonomy_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version 2.5;

use File::Spec::Functions qw(catdir);

sub default_options {
    my ($self) = @_;
    return {
        %{$self->SUPER::default_options},

        user             => $ENV{'USER'},
        email            => $ENV{'USER'} . '@ebi.ac.uk',

        pipeline_name    => 'ncbi_taxonomy_' . $self->o('ensembl_release'),

        taxdump_loc      => 'ftp://ftp.ncbi.nih.gov/pub/taxonomy',
        taxdump_file     => 'taxdump.tar.gz',
        scratch_dir      => catdir('/hps/scratch/flicek/ensembl', $self->o('user'), $self->o('pipeline_name')),

        base_dir         => $ENV{'BASE_DIR'},

    copy_service_uri => $ENV{'DBCOPY_API_URI'},
    src_host         => undef,
    tgt_host         => undef,
    tgt_db_name      => undef,
    copy_to_tgt_host => 0,
    payload          =>
      '{'.
        '"src_host": "'.$self->o('pipeline_db', '-host').':'.$self->o('pipeline_db', '-port').'", '.
        '"src_incl_db": "'.$self->o('user').'_'.$self->o('pipeline_name').'", '.
        '"src_incl_tables": "ncbi_taxa_name,ncbi_taxa_node", '.
        '"tgt_host": "'.$self->o('tgt_host').'", '.
        '"tgt_db_name": "'.$self->o('tgt_db_name').'", '.
        '"user": "'.$self->o('user').'"'.
      '}',
    metadata_host_uri => undef,
    metadata_db_name  => undef,
  };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
        @{$self->SUPER::pipeline_create_commands},
        $self->db_cmd() . ' < ' . $self->o('base_dir') . '/ensembl-taxonomy/sql/table.sql',
        'mkdir -p ' . $self->o('scratch_dir'),
    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},
        taxdump_loc  => $self->o('taxdump_loc'),
        taxdump_file => $self->o('taxdump_file'),
        scratch_dir  => $self->o('scratch_dir'),
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
        'default' => { LSF => '-q production' },
        '16GB'    => { LSF => '-q production -M 16000' },
    };
}

sub pipeline_analyses {
  my ($self) = @_;
  return [
    {
      -logic_name => 'download_tarball',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -input_ids  => [ {} ],
      -parameters => {
                       cmd => 'curl #taxdump_loc#/#taxdump_file# > #scratch_dir#/#taxdump_file#',
                     },
      -flow_into  => { 1 => [ 'untar'], },
    },
    {
      -logic_name => 'untar',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => 'cd #scratch_dir# ; tar -xzf #scratch_dir#/#taxdump_file#',
                     },
      -flow_into  => {
                       '1->A' => [ 'load_nodes' ],
                       'A->1' => [ 'build_left_right_indices' ],
                     },
    },
    {
      -logic_name => 'load_nodes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name    => '16GB',
      -parameters => {
                       inputfile => '#scratch_dir#/nodes.dmp',
                       delimiter => "\t\Q|\E\t?",
                     },
      -flow_into  => {
                       1 => [ 'zero_parent_id', 'uniq_names' ],
                       2 => { '?table_name=ncbi_taxa_node' =>
                              { 'taxon_id' => '#_0#',
                                'parent_id' => '#_1#',
                                'rank' => '#_2#',
                                'genbank_hidden_flag' => '#_10#',
                                'left_index' => 1,
                                'right_index' => 1,
                                'root_id' => 1
                              }
                            },
                      },
    },
    {
      -logic_name => 'zero_parent_id',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
                       sql => "update ncbi_taxa_node set parent_id=0 where parent_id=taxon_id",
                     },
    },
    # This analysis requires the names to be loaded (to find the "root" node)
    {
      -logic_name => 'build_left_right_indices',
      -module     => 'Bio::EnsEMBL::Taxonomy::RunnableDB::AddLeftRightIndexes',
      -flow_into  => { 1 => [ 'add_import_date' ], },
    },
    {
      -logic_name => 'uniq_names',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => 'cut -f1-4,7- #scratch_dir#/names.dmp | uniq > #scratch_dir#/names.uniq.dmp',
                     },
      -flow_into  => { 1 => [ 'load_names' ], },
    },
    {
      -logic_name => 'load_names',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name    => '16GB',
      -parameters => {
                       inputfile => '#scratch_dir#/names.uniq.dmp',
                       delimiter => "\t\Q|\E\t?",
                     },
      -flow_into  => {
                       1 => [ 'load_merged_names' ],
                       2 => { '?table_name=ncbi_taxa_name' =>
                              { 'taxon_id' => '#_0#',
                                'name' => '#_1#',
                                'name_class' => '#_2#'
                              }
                            },
                     },
    },
    {
      -logic_name => 'load_merged_names',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
                       inputfile => '#scratch_dir#/merged.dmp',
                       delimiter => "\t\Q|\E\t?",
                     },
      -flow_into  => {
                       1 => [ 'web_name_patches' ],
                       2 => { '?table_name=ncbi_taxa_name' =>
                              { 'name' => '#_0#',
                                'taxon_id' => '#_1#',
                                'name_class' => 'merged_taxon_id'
                              }
                            },
                      },
    },
    {
      -logic_name    => 'web_name_patches',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters    => {
                          input_file => $self->o('base_dir').'/ensembl-taxonomy/sql/web_name_patches.sql',
                        },
    },
    {
      -logic_name => 'add_import_date',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
                       sql => 'INSERT INTO ncbi_taxa_name (taxon_id, name_class, name) SELECT taxon_id, "import date", CURRENT_TIMESTAMP FROM ncbi_taxa_node WHERE parent_id=0 GROUP BY taxon_id',
                     },
      -flow_into  => { 1 => [ 'cleanup' ], },
    },
    {
      -logic_name => 'cleanup',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => 'rm -rf #scratch_dir#',
                     },
      -flow_into  => { 1 => [ 'post_load_checks' ], },
    },
    {
      -logic_name => 'post_load_checks',
      -module     => 'Bio::EnsEMBL::Taxonomy::RunnableDB::PostLoadChecks',
      -parameters => {
                       tgt_host => $self->o('tgt_host'),
                       copy_to_tgt_host => $self->o('copy_to_tgt_host'),
                     },
      -flow_into  => { 1 => WHEN(
                                    'defined #tgt_host# && #copy_to_tgt_host#' => [ 'copy_database' ],
                                ),
                     },
    },
    {
      -logic_name => 'copy_database',
      -module     => 'ensembl.production.hive.ProductionDBCopy',
      -language   => 'python3',
      -rc_name    => 'default',
      -parameters => {
                       endpoint => $self->o('copy_service_uri'),
                       payload  => $self->o('payload'),
                       method   => 'post',
                    },
      -flow_into  => { 1 => WHEN(
                            'defined #metadata_host_uri# && #metadata_db_name#'  => ['update_taxonomy_node_table_in_metadata_db']
                        ),
      },
    },
    {
      -logic_name => 'update_taxonomy_node_table_in_metadata_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
                       db_conn => $self->o('metadata_host_uri').$self->o('metadata_db_name'),
                       sql => qq{ INSERT INTO } . $self->o('metadata_db_name'). qq{.ncbi_taxa_node (taxon_id, parent_id,
                                  rank, genbank_hidden_flag, left_index, right_index, root_id)
                                  SELECT taxon_id, parent_id, rank, genbank_hidden_flag, left_index, right_index, root_id
                                  from } . $self->o('tgt_db_name') .
                                  qq{.ncbi_taxa_node t2 order by taxon_id
                                  on duplicate key update taxon_id = t2.taxon_id, parent_id = t2.parent_id, rank = t2.rank,
                                  genbank_hidden_flag=t2.genbank_hidden_flag, left_index=t2.left_index, right_index=t2.right_index,
                                  root_id=t2.root_id;}
                      },
      -flow_into  => { 1 => 'update_taxonomy_name_table_in_metadata_db' },
    },
    {
      -logic_name => 'update_taxonomy_name_table_in_metadata_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
                       db_conn => $self->o('metadata_host_uri').$self->o('metadata_db_name'),
                       sql => qq{ INSERT INTO } . $self->o('metadata_db_name') . qq{.ncbi_taxa_name (taxon_id, name, name_class)
                            SELECT taxon_id, name, name_class from } . $self->o('tgt_db_name') . qq{.ncbi_taxa_name t2 order by taxon_id
                            on duplicate key update taxon_id = t2.taxon_id, name = t2.name, name_class = t2.name_class;
                       }
                      },
    },
  ];
}

1;
