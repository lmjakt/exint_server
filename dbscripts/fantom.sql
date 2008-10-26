/* create some tables to keep fantom data in. Should probably include the sequence as well.. not too tricky..

  this might as well include the sequence as we have one sequence for each transcript.. 
*/

create table fantom_transcripts (
	transcript int unique,
	id text,
	length int,
	sequence text
);

/* a table for the annotation fields .. */
create table fantom_fields (
	field int unique,
	description text
);

/* create a table for the annotation */
create table fantom_annotation (
	transcript int references fantom_transcripts(transcript),
	field int references fantom_fields(field),
	annotation text,
	dbref text,
	evidence text
);

/* create a table for the mapping of the transcripts onto the mouse genome */
create table fantom_assemblies (
	assembly int unique,
	transcript int references fantom_transcripts(transcript),
	chromosome text,
	strand int,
	score float,
	rank int
);

create table fantom_matches (
	assembly int references fantom_assemblies (assembly),
	cstart int,
	cend int,
	fstart int,
	fend int,
	percent float
);

