# frozen_string_literal: true

module Bio
  class VCF
    # Rewrites GT and Number=G fields for one selected ALT allele.
    class GenotypeRewriter
      def split_gt(value, selected_index)
        return value if value.nil? || value == '.'

        tokens = value.scan(/[0-9.]+|[\/|]/)
        alleles = even_values(tokens)
        separators = odd_values(tokens)
        mapped = alleles.collect { |allele| map_allele(allele, selected_index) }
        join_genotype(mapped, separators)
      end

      def split_values(values, allele_count, selected_index, ploidy)
        return values unless ploidy

        combinations(allele_count, ploidy).each_with_index.collect do |combination,
                                                                   index|
          next unless retained?(combination, selected_index)

          values.fetch(index)
        end.compact
      end

      def split_pl(values, allele_count, selected_index, ploidy)
        selected = split_values(values, allele_count, selected_index, ploidy)
        normalize_pl(selected)
      end

      def normalize_selected_pl(values)
        normalize_pl(values).join(',')
      end

      def ploidy(value)
        return nil unless value && value != '.'

        value.scan(/[0-9.]+/).length
      end

      def infer_ploidy(value_count, allele_count)
        ploidy = 1
        loop do
          count = genotype_count(allele_count, ploidy)
          return ploidy if count == value_count
          return nil if count > value_count

          ploidy += 1
        end
      end

      private

      def map_allele(allele, selected_index)
        return '.' if allele == '.'

        numeric = allele.to_i
        numeric.zero? ? '0' : (numeric == selected_index + 1 ? '1' : '.')
      end

      def join_genotype(alleles, separators)
        alleles.each_with_index.collect do |allele, index|
          allele + (separators[index] || String.new)
        end.join
      end

      def retained?(combination, selected_index)
        combination.all? { |allele| allele.zero? || allele == selected_index + 1 }
      end

      def combinations(allele_count, ploidy)
        (0...allele_count).to_a.repeated_combination(ploidy).to_a
      end

      def genotype_count(allele_count, ploidy)
        combinations(allele_count, ploidy).length
      end

      def normalize_pl(values)
        numeric_values = values.reject { |value| value == '.' }.collect(&:to_i)
        minimum = numeric_values.min
        values.collect do |value|
          value == '.' ? value : value.to_i - minimum
        end
      end

      def even_values(values)
        values.values_at(*values.each_index.select(&:even?))
      end

      def odd_values(values)
        values.values_at(*values.each_index.select(&:odd?))
      end
    end
  end
end
