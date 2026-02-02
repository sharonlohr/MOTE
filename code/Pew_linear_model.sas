PROC IMPORT OUT= WORK.mercer 
            DATAFILE= "C:\mercer.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

PROC GLM DATA=MERCER;
  CLASS TYPE (ref="prob") TOPIC (ref="Other") SAMPID;
  MODEL emp_biassq_vi = TYPE topic type*topic/solution;
  LSMEANS TYPE/DIFF;
  LSMEANS TYPE*TOPIC/DIFF;
  title 'Linear model, no random effect';
run;


PROC MIXED DATA=MERCER;
  CLASS TYPE (ref="prob") TOPIC (ref="Other") SAMPID;
  MODEL EMP_BIASSQ_VI = TYPE TOPIC TYPE*TOPIC/SOLUTION;
  RANDOM SAMPID(TYPE);
  LSMEANS TYPE/DIFF;
  LSMEANS TYPE*TOPIC/DIFF;
  title 'Linear mixed model with main effects and interactions';
run;


PROC MIXED DATA=MERCER;
  CLASS TYPE (ref="prob") TOPIC (ref="Other") SAMPID;
  MODEL emp_biassq_vi =  type*topic/solution;
  RANDOM sampid(TYPE);
  LSMEANS TYPE*TOPIC/DIFF;
  TITLE 'Linear mixed model done as one-way';
run;

