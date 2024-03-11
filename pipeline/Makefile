all: tune_wflows merge_results summarize_results

splendid: train_eval merge_eval iv_summary

sequential: seq_tune seq_merge seq_summarize

retrain: retrain_tune retrain_merge retrain_summarize

gene_opt: gene_tune gene_seq_train_eval gene_seq_merge_eval

tune_wflows:
	./1-tune_wflows.sh $(filter-out $@,$(MAKECMDGOALS))

merge_results:
	./2-merge_results.sh $(filter-out $@,$(MAKECMDGOALS))

summarize_results:
	./3-summarize_results.sh $(filter-out $@,$(MAKECMDGOALS))

# Splendid pipeline
# Train models and return evaluations
train_eval:
	./splendid/1-train_eval.sh $(filter-out $@,$(MAKECMDGOALS))

# Merge evaluations across reps and algorithms
merge_eval:
	./splendid/2-merge_eval.sh $(filter-out $@,$(MAKECMDGOALS))

# Internal validation summary and combined across datasets
iv_summary:
	./splendid/3-iv_summary.sh $(filter-out $@,$(MAKECMDGOALS))

# Sequential classifier
seq_tune:
	./sequential/1-tune_wflows.sh $(filter-out $@,$(MAKECMDGOALS))

seq_merge:
	./sequential/2-merge_results.sh $(filter-out $@,$(MAKECMDGOALS))

seq_summarize:
	./sequential/3-summarize_results.sh $(filter-out $@,$(MAKECMDGOALS))

# Retrain for sequential and two-step classifiers
retrain_tune:
	./retrain/1-tune_wflows.sh $(filter-out $@,$(MAKECMDGOALS))

retrain_merge:
	./retrain/2-merge_results.sh $(filter-out $@,$(MAKECMDGOALS))

retrain_summarize:
	./retrain/3-summarize_results.sh $(filter-out $@,$(MAKECMDGOALS))

# Gene optimization
gene_tune:
	./gene_opt/1a-tune_wflows.sh $(filter-out $@,$(MAKECMDGOALS))

gene_merge_eval:
	./gene_opt/2a-gene_merge_eval.sh $(filter-out $@,$(MAKECMDGOALS))

gene_seq_train_eval:
	./gene_opt/1b-gene_seq_train_eval.sh $(filter-out $@,$(MAKECMDGOALS))

gene_seq_merge_eval:
	./gene_opt/2b-gene_seq_merge_eval.sh $(filter-out $@,$(MAKECMDGOALS))

# Clean target
%:
		@:
clean:
	./assets/clean.sh
