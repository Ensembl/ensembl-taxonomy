-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2024] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

/**
@table ncbi_taxa_node
@desc This table contains all taxa used in this database, which mirror the data and tree structure from NCBI Taxonomy database
@colour   #3CB371

@example    This examples shows how to get the lineage for Homo sapiens:
    @sql    SELECT n2.taxon_id, n2.parent_id, na.name, n2.rank, n2.left_index, n2.right_index FROM ncbi_taxa_node n1 JOIN (ncbi_taxa_node n2 LEFT JOIN ncbi_taxa_name na ON n2.taxon_id = na.taxon_id AND na.name_class = "scientific name")  ON n2.left_index <= n1.left_index AND n2.right_index >= n1.right_index WHERE n1.taxon_id = 9606 ORDER BY left_index;

@column taxon_id                The NCBI Taxonomy ID
@column parent_id               The parent taxonomy ID for this node (refers to ncbi_taxa_node.taxon_id)
@column rank                    E.g. kingdom, family, genus, etc.
@column genbank_hidden_flag     Boolean value which defines whether this rank is used or not in the abbreviated lineage
@column left_index              Sub-set left index. All sub-nodes have left_index and right_index values larger than this left_index
@column right_index             Sub-set right index. All sub-nodes have left_index and right_index values smaller than this right_index
@column root_id                 The root taxonomy ID for this node (refers to ncbi_taxa_node.taxon_id)

@see ncbi_taxa_name
*/

CREATE TABLE ncbi_taxa_node (
  taxon_id                        int(10) unsigned NOT NULL,
  parent_id                       int(10) unsigned NOT NULL,

  rank                            char(32) default '' NOT NULL,
  genbank_hidden_flag             tinyint(1) default 0 NOT NULL,

  left_index                      int(10) DEFAULT 0 NOT NULL,
  right_index                     int(10) DEFAULT 0 NOT NULL,
  root_id                         int(10) default 1 NOT NULL,

  PRIMARY KEY (taxon_id),
  KEY (parent_id),
  KEY (rank),
  KEY (left_index),
  KEY (right_index)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

/**
@table ncbi_taxa_name
@desc This table contains different names, aliases and meta data for the taxa used in Ensembl.
@colour   #3CB371

@example    Here is an example on how to get the taxonomic ID for a species:
    @sql                          SELECT * FROM ncbi_taxa_name WHERE name_class = "scientific name" AND name = "Homo sapiens";

@column taxon_id              External reference to taxon_id in @link ncbi_taxa_node
@column name                  Information assigned to this taxon_id
@column name_class            Type of information. e.g. common name, genbank_synonym, scientif name, etc.

@see ncbi_taxa_node
*/


CREATE TABLE ncbi_taxa_name (
  taxon_id                    int(10) unsigned NOT NULL,
  name                        varchar(500) NOT NULL,
  name_class                  varchar(50) NOT NULL,

  FOREIGN KEY (taxon_id) REFERENCES ncbi_taxa_node(taxon_id),

	-- NO PK

  KEY (taxon_id),
  KEY (name),
  KEY (name_class)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


