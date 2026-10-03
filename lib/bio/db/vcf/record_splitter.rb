# frozen_string_literal: true

module Bio
  class VCF
    # Splits one multiallelic record into biallelic records.
    class RecordSplitter
      def initialize
        @rewriter = FieldRewriter.new
      end

      def split(record, selected_index)
        result = record.dup
        result.info = split_info(record, selected_index)
        result.samples = split_samples(record, selected_index)
        result.alt = [record.alt[selected_index]]
        result.normalization_status = nil
        result
      end

      private

      def split_info(record, selected_index)
        record.info.each_with_object({}) do |(key, value), hash|
          hash[key] = @rewriter.split_value(value, number(record, :info, key),
            selected_index, record.alt.length + 1, nil)
        end
      end

      def split_samples(record, selected_index)
        record.samples.each_with_object({}) do |(name, sample), hash|
          hash[name] = split_sample(record, sample, selected_index)
        end
      end

      def split_sample(record, sample, selected_index)
        allele_count = record.alt.length + 1
        ploidy = @rewriter.ploidy_from_gt(sample['GT'])
        sample.each_with_object({}) do |(key, value), result|
          result[key] = split_field(record, key, value, selected_index,
                                    allele_count, ploidy)
        end
      end

      def split_field(record, key, value, selected_index, allele_count, ploidy)
        return @rewriter.split_gt(value, selected_index) if key == 'GT'
        return @rewriter.split_pl(value, allele_count, selected_index,
                                 ploidy) if key == 'PL'

        @rewriter.split_value(value, number(record, :format, key),
          selected_index, allele_count, ploidy)
      end

      def number(record, kind, key)
        header = record.header
        definitions = header && (kind == :info ? header.info : header.format)
        definition = definitions && definitions[key]
        definition ? definition['Number'] : builtin_number(kind, key)
      end

      def builtin_number(kind, key)
        if kind == :info
          FieldRewriter::INFO_NUMBERS[key]
        else
          FieldRewriter::FORMAT_NUMBERS[key]
        end
      end
    end
  end
end
