# frozen_string_literal: true

module Bio
  class VCF
    # Parses one VCF data line into a Record.
    module RecordParser
      def parse(line)
        columns = line.chomp.split("\t", -1)
        @chrom, pos_string, @id, @ref, alt_string, @qual, @filter,
          info_string, format_string = columns
        @pos = pos_string.to_i
        @alt = alt_string == '.' ? [] : alt_string.split(',')
        parse_info(info_string)
        parse_samples(columns, format_string)
      end

      private

      def parse_info(info_string)
        return if info_string.nil? || info_string == '.'

        info_string.split(';').each do |item|
          key, value = item.split('=', 2)
          @info[key] = value.nil? ? true : value
        end
      end

      def parse_samples(columns, format_string)
        return if format_string.nil? || format_string == '.'

        @format = format_string.split(':')
        names = header ? header.sample_names : []
        columns[9..-1].to_a.each_with_index do |sample_string, index|
          name = names[index] || "SAMPLE#{index + 1}"
          values = sample_string.split(':', -1)
          @samples[name] = Hash[format.zip(values)]
        end
      end
    end
  end
end
