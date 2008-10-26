/* Add some tables that hold annotation from affymetrix in a
   more organised and flexible way than before. Basically add
   a two tables, one containing fields and the other containing
   annotation, as for the ensemb annoation
*/

-- first a sequence for the fields..
create sequence affymetrix_annotation_fields_seq;
create table affymetrix_annotation_fields (
affy_field int default nextval('affymetrix_annotation_fields_seq') primary key,
field_name text);

-- and then the actual annotation
create sequence affymetrix_annotation_seq;
create table affymetrix_annotation (
annotation_index int default nextval('affymetrix_annotation_seq') primary key,
index int references p_sets (index),
affy_field int references affymetrix_annotation_fields (affy_field),
annotation text);

-- which should be ok. maybe
