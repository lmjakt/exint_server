--- add a table that stores the individual affymetrix probes and sequences
create table probes (
probe serial primary key,
p_set int references p_sets(index),
x_pos int,
y_pos int,
probe_pos int,
seq text,
strand text);

