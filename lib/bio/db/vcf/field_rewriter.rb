# frozen_string_literal: true

module Bio
  class VCF
    # Rewrites Number=A, Number=R, and Number=G fields after ALT splitting.
    class FieldRewriter
      INFO_NUMBERS = { 'AC' => 'A', 'AF' => 'A' }.freeze
      FORMAT_NUMBERS = {
        'GT' => '1', 'AD' => 'R', 'ADF' => 'R', 'ADR' => 'R',
        'GL' => 'G', 'PL' => 'G', 'PP' => 'G'
      }.freeze

      def initialize
        @genotype_rewriter = GenotypeRewriter.new
      end

      def split_value(value, number, selected_index, allele_count, ploidy)
        return value unless value.is_a?(String) && value != '.' && number

        case number
        when 'A' then select_values(value, [selected_index])
        when 'R' then select_values(value, [0, selected_index + 1])
        when 'G' then split_genotype(value, allele_count, selected_index,
                                     ploidy)
        else value
        end
      rescue IndexError
        value
      end

      def split_gt(value, selected_index)
        @genotype_rewriter.split_gt(value, selected_index)
      end

      def split_pl(value, allele_count, selected_index, ploidy)
        return value unless value.is_a?(String) && value != '.' && ploidy

        values = value.split(',', -1)
        selected = @genotype_rewriter.split_values(values, allele_count,
                                                   selected_index, ploidy)
        return value if selected.empty?

        @genotype_rewriter.normalize_selected_pl(selected)
      rescue IndexError
        value
      end

      def ploidy_from_gt(value)
        @genotype_rewriter.ploidy(value)
      end

      private

      def select_values(value, indexes)
        value.split(',', -1).values_at(*indexes).join(',')
      end

      def split_genotype(value, allele_count, selected_index, ploidy)
        values = value.split(',', -1)
        ploidy ||= @genotype_rewriter.infer_ploidy(values.length, allele_count)
        return value unless ploidy

        selected = @genotype_rewriter.split_values(values, allele_count,
                                                   selected_index, ploidy)
        selected.empty? ? value : selected.join(',')
      end
    end
  end
end
