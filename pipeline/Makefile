all: train_eval merge_eval iv_summary

sequential: seq_train_eval seq_merge_eval seq_iv_summary

retrain: retrain_eval retrain_merge retrain_summary

gene_opt: gene_train_eval gene_seq_train_eval

# Train models and return evaluations
train_eval:
	./1-train_eval.sh $(filter-out $@,$(MAKECMDGOALS))

# Merge evaluations across reps and algorithms
merge_eval:
	./2-merge_eval.sh $(filter-out $@,$(MAKECMDGOALS))

# Internal validation summary and combined across datasets
iv_summary:
	./3-iv_summary.sh $(filter-out $@,$(MAKECMDGOALS))

# Sequential classifier
seq_train_eval:
	./sequential/1-seq_train_eval.sh $(filter-out $@,$(MAKECMDGOALS))

seq_merge_eval:
	./sequential/2-seq_merge_eval.sh $(filter-out $@,$(MAKECMDGOALS))

seq_iv_summary:
	./sequential/3-seq_iv_summary.sh $(filter-out $@,$(MAKECMDGOALS))

# Retrain for sequential and two-step classifiers
retrain_eval:
	./retrain/1-retrain_eval.sh $(filter-out $@,$(MAKECMDGOALS))

retrain_merge:
	./retrain/2-retrain_merge.sh $(filter-out $@,$(MAKECMDGOALS))

retrain_summary:
	./retrain/3-retrain_summary.sh $(filter-out $@,$(MAKECMDGOALS))

# Gene optimization
gene_train_eval:
	./gene_opt/1a-gene_train_eval.sh $(filter-out $@,$(MAKECMDGOALS))

gene_seq_train_eval:
	./gene_opt/1b-gene_seq_train_eval.sh $(filter-out $@,$(MAKECMDGOALS))

# Clean target
%:
		@:
clean:
	./assets/clean.sh
