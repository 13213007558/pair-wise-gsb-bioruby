# frozen_string_literal: true

module Bio
  class VCF
    # Serializes a Record back to VCF text.
    module RecordSerializer
      def to_s
        (fixed_columns + sample_columns).join("\t") + "\n"
      end

      private

      def fixed_columns
        [
          chrom, pos.to_s, value_or_dot(id), value_or_dot(ref),
          alt_to_s, value_or_dot(qual), value_or_dot(filter), info_to_s
        ]
      end

      def sample_columns
        return [] unless format.any?

        names = header ? header.sample_names : samples.keys
        [format.join(':')] + names.collect { |name| sample_to_s(name) }
      end

      def value_or_dot(value)
        value.nil? || value.to_s.empty? ? '.' : value.to_s
      end

      def alt_to_s
        alt.empty? ? '.' : alt.join(',')
      end

      def info_to_s
        return '.' if info.empty?

        info.collect do |key, value|
          value == true ? key : "#{key}=#{value}"
        end.join(';')
      end

      def sample_to_s(name)
        sample = samples[name]
        return '.' unless sample

        format.collect { |key| sample[key] || '.' }.join(':')
      end
    end
  end
end
