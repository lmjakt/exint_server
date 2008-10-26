/*
	Note that you'll need to change the n value to something reasonable for your
	system. max(index) + 1 is a reasonable number to set it to, but I don't know
	how to write one expression that does that in SQL
*/
create sequence cel_files_index_seq;
select setval('cel_files_index_seq', n);
alter table cel_files alter column index set default nextval('cel_files_index_seq');
alter table cel_files add unique(index);

/*     and the same for the experiments table.. */
create sequence experiments_index_seq;
select setval('experiments_index_seq', n);
alter table experiments alter column index set default nextval('experiments_index_seq');
alter table experiments add unique(index);

