
module Bio
  module Phylo
    # Log-likelihood calculation for nucleotide alignments under JC69.
    #
    # Leaf names must exactly match the keys of a hash-like alignment such as
    # a Hash, Bio::Alignment::SequenceHash, or Bio::Alignment::OriginalAlignment.
    # IUPAC ambiguity characters denote the set of allowed states. Gap (
    # <tt>-</tt>, <tt>.</tt>) and missing-data (<tt>?</tt>) characters are
    # treated as unobserved rather than as a fifth state. A site at which all
    # leaves are unobserved therefore has log-likelihood 0.
    #
    # Calculations are performed in log space. Multifurcating trees and
    # zero-length branches are supported, and long trees do not underflow to
    # negative infinity from multiplying small probabilities.
    class JC69
      IUPAC_STATES = {
        'A' => [0].freeze,
        'C' => [1].freeze,
        'G' => [2].freeze,
        'T' => [3].freeze,
        'U' => [3].freeze,
        'M' => [0, 1].freeze,
        'R' => [0, 2].freeze,
        'W' => [0, 3].freeze,
        'S' => [1, 2].freeze,
        'Y' => [1, 3].freeze,
        'K' => [2, 3].freeze,
        'V' => [0, 1, 2].freeze,
        'H' => [0, 1, 3].freeze,
        'D' => [0, 2, 3].freeze,
        'B' => [1, 2, 3].freeze,
        'N' => [0, 1, 2, 3].freeze
      }.freeze

      MISSING_CHARACTERS = ['-', '.', '?'].freeze
      DEFAULT_RATE_CATEGORIES = [[1.0, 1.0].freeze].freeze
      LOG_ZERO = -Float::INFINITY
      ROOT_LOG_PROBABILITY = -Math.log(4.0)

      # Creates and runs a JC69 likelihood calculation.
      #
      # Each rate category is represented as [rate, weight]. Non-negative
      # weights do not have to sum to one; positive weights are normalized.
      def self.calculate(tree, alignment, rate_categories = DEFAULT_RATE_CATEGORIES)
        new(tree, alignment, rate_categories)
      end

      # Returns the total log-likelihood.
      def self.log_likelihood(tree, alignment, rate_categories = DEFAULT_RATE_CATEGORIES)
        new(tree, alignment, rate_categories).log_likelihood
      end

      # Returns the log-likelihood of every alignment column.
      def self.site_log_likelihoods(tree, alignment,
                                    rate_categories = DEFAULT_RATE_CATEGORIES)
        new(tree, alignment, rate_categories).site_log_likelihoods
      end

      attr_reader :site_log_likelihoods

      def initialize(tree, alignment, rate_categories = DEFAULT_RATE_CATEGORIES)
        @tree = tree
        @alignment = alignment
        @rates, @log_weights = prepare_rate_categories(rate_categories)

        validate_inputs
        @root, @children, @parents, @postorder = rooted_tree
        @sequences = extract_sequences
        validate_alignment_labels
        @alignment_length = sequence_length
        @branch_log_probabilities = build_branch_log_probabilities
        @site_log_likelihoods = calculate_site_log_likelihoods
        @log_likelihood = @site_log_likelihoods.sum
      end

      # Total log-likelihood, summed over alignment columns and rate categories.
      def log_likelihood
        @log_likelihood
      end

      private

      def prepare_rate_categories(categories)
        categories = categories || DEFAULT_RATE_CATEGORIES
        unless categories.respond_to?(:each) && !categories.is_a?(Hash)
          raise ArgumentError, 'rate categories must be an array of [rate, weight]'
        end
        raise ArgumentError, 'at least one rate category is required' if categories.empty?

        prepared = categories.each_with_index.map do |category, index|
          rate, weight = rate_and_weight(category, index)
          validate_rate(rate, index)
          validate_weight(weight, index)
          [rate.to_f, weight.to_f]
        end
        prepared.select! { |_rate, weight| weight > 0 }
        if prepared.empty?
          raise ArgumentError, 'rate category weights must have a positive sum'
        end

        weight_sum = prepared.sum { |_rate, weight| weight }
        log_weights = prepared.map do |_rate, weight|
          Math.log(weight / weight_sum)
        end
        [prepared.map(&:first), log_weights]
      end

      def rate_and_weight(category, index)
        if category.is_a?(Hash)
          [category.fetch(:rate), category.fetch(:weight)]
        elsif category.respond_to?(:to_ary)
          values = category.to_ary
          unless values.length == 2
            raise ArgumentError, "rate category #{index} must contain a rate and weight"
          end
          values
        elsif category.respond_to?(:rate) && category.respond_to?(:weight)
          [category.rate, category.weight]
        else
          raise ArgumentError, "rate category #{index} must be [rate, weight]"
        end
      rescue KeyError => error
        raise ArgumentError, "rate category #{index} is missing #{error.key}"
      end

      def validate_rate(rate, index)
        unless rate.is_a?(Numeric) && rate.finite? && rate >= 0
          raise ArgumentError, "rate #{index} must be a finite non-negative number"
        end
      end

      def validate_weight(weight, index)
        unless weight.is_a?(Numeric) && weight.finite? && weight >= 0
          raise ArgumentError, "weight #{index} must be a finite non-negative number"
        end
      end

      def validate_inputs
        unless @tree.is_a?(Bio::Tree)
          raise ArgumentError, 'tree must be a Bio::Tree'
        end
        unless @alignment.respond_to?(:each_pair)
          raise ArgumentError, 'alignment must provide labelled sequences via each_pair'
        end
        raise ArgumentError, 'tree must not be empty' if @tree.nodes.empty?
      end

      def rooted_tree
        nodes = @tree.nodes
        root = nodes.find { |node| @tree.out_degree(node) != 1 } || nodes.first

        parents = { root => nil }
        visited = {}
        preorder = []
        stack = [root]

        until stack.empty?
          node = stack.pop
          raise ArgumentError, 'tree must be connected and acyclic' if visited.key?(node)
          visited[node] = true
          preorder << node

          @tree.adjacent_nodes(node).each do |neighbor|
            next if neighbor == parents[node]
            if parents.key?(neighbor)
              raise ArgumentError, 'tree must be connected and acyclic'
            end
            parents[neighbor] = node
            stack << neighbor
          end
        end

        unless visited.length == nodes.length && @tree.edges.length == nodes.length - 1
          raise ArgumentError, 'tree must be connected and acyclic'
        end

        children = Hash.new { |hash, key| hash[key] = [] }
        parents.each do |node, parent|
          children[parent] << node unless parent.nil?
        end
        parents.default = nil
        [root, children, parents, preorder.reverse]
      end

      def extract_sequences
        sequences = {}
        @alignment.each_pair do |label, sequence|
          name = String(label)
          if sequences.key?(name)
            raise ArgumentError, "duplicate alignment label #{name.inspect}"
          end
          unless sequence.respond_to?(:to_str)
            raise ArgumentError, "alignment sequence #{name.inspect} must be a string"
          end
          sequences[name] = sequence.to_str
        end
        sequences
      end

      def validate_alignment_labels
        @terminal_nodes = @tree.nodes.select do |node|
          @tree.out_degree(node) <= 1
        end
        terminal_labels = @terminal_nodes.map do |node|
          name = @tree.get_node_name(node)
          name.nil? ? nil : String(name)
        end

        if terminal_labels.any?(&:nil?) || terminal_labels.any?(&:empty?)
          raise ArgumentError, 'every leaf must have a non-empty label'
        end
        if terminal_labels.uniq.length != terminal_labels.length
          raise ArgumentError, 'leaf labels must be unique'
        end

        missing = terminal_labels - @sequences.keys
        extra = @sequences.keys - terminal_labels
        unless missing.empty? && extra.empty?
          message = +'alignment labels do not match leaf labels'
          message << "; missing: #{missing.sort.inspect}" unless missing.empty?
          message << "; extra: #{extra.sort.inspect}" unless extra.empty?
          raise ArgumentError, message
        end
      end

      def sequence_length
        lengths = @sequences.values.map(&:length)
        if lengths.uniq.length > 1
          raise ArgumentError, 'alignment sequences must have the same length'
        end
        lengths.first || 0
      end

      def build_branch_log_probabilities
        @postorder.each_with_object({}) do |node, probabilities|
          next if node == @root

          parent = @parents[node]
          edge = @tree.get_edge(parent, node)
          length = @tree.get_edge_distance(edge)
          unless length.is_a?(Numeric) && length.finite? && length >= 0
            raise ArgumentError, 'every branch length must be a finite non-negative number'
          end

          probabilities[node] = @rates.map do |rate|
            transition_log_probabilities(rate * length)
          end
        end
      end

      def transition_log_probabilities(scaled_length)
        scaled_length = scaled_length.to_f
        survival = if scaled_length.infinite?
                     0.0
                   else
                     Math.exp(-4.0 * scaled_length / 3.0)
                   end
        same = 0.25 + 0.75 * survival
        different = 0.25 - 0.25 * survival
        different = 0.0 if different.negative?
        [Math.log(same), different.zero? ? LOG_ZERO : Math.log(different)]
      end

      def calculate_site_log_likelihoods
        Array.new(@alignment_length) do |position|
          next 0.0 if uninformative_site?(position)

          category_log_likelihoods = @rates.each_with_index.map do |_rate, rate_index|
            category_log_likelihood(position, rate_index)
          end
          log_sum_exp(category_log_likelihoods.zip(@log_weights).map do |likelihood, log_weight|
            likelihood + log_weight
          end)
        end
      end

      def uninformative_site?(position)
        @terminal_nodes.all? do |node|
          states = character_states(node, position)
          states.nil? || states.length == 4
        end
      end

      def category_log_likelihood(position, rate_index)
        partial_log_likelihoods = {}

        @postorder.each do |node|
          children = @children[node]
          partial_log_likelihoods[node] = if children.empty?
                                            leaf_vector(node, position)
                                          else
                                            internal_vector(children,
                                                            partial_log_likelihoods, rate_index)
                                          end
        end

        root_vector = partial_log_likelihoods[@root]
        if @tree.out_degree(@root) <= 1
          own_vector = leaf_vector(@root, position)
          root_vector = root_vector.each_with_index.map do |likelihood, state|
            likelihood + own_vector[state]
          end
        end

        root_terms = root_vector.each_with_index.map do |likelihood, _state|
          likelihood + ROOT_LOG_PROBABILITY
        end
        log_sum_exp(root_terms)
      end

      def leaf_vector(node, position)
        states = character_states(node, position)
        Array.new(4) do |state|
          states.nil? || states.include?(state) ? 0.0 : LOG_ZERO
        end
      end

      def internal_vector(children, partial_log_likelihoods, rate_index)
        Array.new(4) do |state|
          children.sum do |child|
            child_likelihoods = partial_log_likelihoods[child]
            terms = child_likelihoods.each_with_index.map do |likelihood, child_state|
              transition = @branch_log_probabilities[child][rate_index]
              probability = state == child_state ? transition.first : transition.last
              likelihood + probability
            end
            log_sum_exp(terms)
          end
        end
      end

      def character_states(node, position)
        label = String(@tree.get_node_name(node))
        character = @sequences[label][position, 1].to_s.upcase
        return nil if missing_character?(character)

        states = IUPAC_STATES[character]
        if states.nil?
          message = "unsupported nucleotide character #{character.inspect} " \
                    "at #{label.inspect}, site #{position}"
          raise ArgumentError, message
        end
        states
      end

      def missing_character?(character)
        MISSING_CHARACTERS.include?(character) ||
          (@alignment.respond_to?(:gap_char) && character == @alignment.gap_char.to_s.upcase) ||
          (@alignment.respond_to?(:missing_char) && character == @alignment.missing_char.to_s.upcase)
      end

      def log_sum_exp(log_values)
        finite_values = log_values.reject { |value| value == LOG_ZERO }
        return LOG_ZERO if finite_values.empty?
        maximum = finite_values.max
        sum = finite_values.sum { |value| Math.exp(value - maximum) }
        maximum + Math.log(sum)
      end
    end
  end
end
