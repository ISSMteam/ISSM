#ifndef PERF_HPP
#define PERF_HPP

#include <chrono>
#include <map>
#include <string>

#include "../../shared/io/Comm/IssmComm.h"
#include "../mpi/issmmpi.h"

struct Perf {

    struct DataEntry {
        std::string marker = "";

        int calls = 0;
        double min_time = 0.0;
        double max_time = 0.0;
        double total_time = 0.0;

        void init(std::string const& m, double t) {
            marker = m;
            calls = 1;
            min_time = t;
            max_time = t;
            total_time = t;
        }

        void add(double t) {
            calls += 1;

            min_time = std::min(min_time, t);
            max_time = std::max(max_time, t);
            total_time += t;
        }

        void add(DataEntry const& o) {
            calls += o.calls;

            min_time = std::min(min_time, o.min_time);
            max_time = std::max(max_time, o.max_time);
            total_time += o.total_time;
        }

        int byteSize() const {
            return sizeof(int) * 2  // calls + marker size
                   + sizeof(double) * 3
                   + (int)(marker.size() + 1); // +1 for null terminator
        }

        char* writeBytes(char* data) const {
            data = write<int>(data, calls);
            data = write<double>(data, min_time);
            data = write<double>(data, max_time);
            data = write<double>(data, total_time);
            data = write<std::string>(data, marker);

            return data;
        }

        template<typename T>
        static char* write(char* data, T const& v) {
            if constexpr (std::is_same_v<T, std::string>) {
                int size = (int)v.size();
                data = write<int>(data, size + 1);
                for(int i = 0; i < size; i += 1) {
                    data = write<char>(data, v[i]);
                }
                data = write<char>(data, 0);

                return data;
            }
            else {
                T* data_cast = reinterpret_cast<T*>(data);
                *data_cast = v;
                return data + sizeof(T);
            }
        }

        char* readBytes(char* data) {
            data = read<int>(data, calls);
            data = read<double>(data, min_time);
            data = read<double>(data, max_time);
            data = read<double>(data, total_time);
            data = read<std::string>(data, marker);

            return data;
        }

        template<typename T>
        static char* read(char* data, T& v) {
            if constexpr (std::is_same_v<T, std::string>) {
                int size = 0;
                data = read<int>(data, size);

                v = data;
                return data + size;
            }
            else {
                T* data_cast = reinterpret_cast<T*>(data);
                v = *data_cast;
                return data + sizeof(T);
            }
        }

        static void writeFileHeader(std::ofstream& out) {
            out << "marker; calls; total_time; avg_time; min_time; max_time\n";
        }

        void writeFile(std::ofstream& out) const {
            double avg_time = total_time / (double)calls;
            bool first = true;
            auto write_func = [&](auto const& v) {
                if(!first) {
                    out << "; ";
                }
                first = false;
                out << v;
            };

            write_func(marker);
            write_func(calls);
            write_func(total_time);
            write_func(avg_time);
            write_func(min_time);
            write_func(max_time);
            out << "\n";
        }
    };


    std::map<std::string, DataEntry> entries;

    void add(std::string const& marker, double time) {
        DataEntry& entry = entries[marker];

        if(0 == entry.calls) {
            entry.init(marker, time);
        } else {
            entry.add(time);
        }
    }

    void write(std::string const& file) {
        int total_ranks = IssmComm::GetSize();

        int my_size = 0;
        for(auto const& cur : entries) {
            my_size += cur.second.byteSize();
        }

        std::vector<int> all_sizes(total_ranks);
        ISSM_MPI_Gather(&my_size, 1, ISSM_MPI_INT, all_sizes.data(), 1, ISSM_MPI_INT, 0, IssmComm::GetComm());

        int total_size = 0;
        std::vector<char> all_data(0);
        std::vector<int> displs(0);
        if(0 == IssmComm::GetRank()) {
            displs.resize(total_ranks);
            for(int i = 0; i < total_ranks; i += 1 ) {
                displs[i] = total_size;
                total_size += all_sizes[i];
            }
            all_data.resize(total_size);
        }

        std::vector<char> local_data(my_size);
        char* pos = local_data.data();
        for(auto const& cur : entries) {
            pos = cur.second.writeBytes(pos);
        }
        assert(my_size == (pos - local_data.data()));

        ISSM_MPI_Gatherv(local_data.data(), my_size, ISSM_MPI_CHAR, all_data.data(), all_sizes.data(), displs.data(), ISSM_MPI_CHAR, 0, IssmComm::GetComm());

        if (0 == IssmComm::GetRank()) {
            std::map<std::string, DataEntry> all_entries;

            pos = all_data.data();
            DataEntry cur_entry = {};
            while(pos < all_data.data() + all_data.size()) {
                pos = cur_entry.readBytes(pos);

                addAll(all_entries, cur_entry);
            }

            std::ofstream out(file);
            DataEntry::writeFileHeader(out);
            for(auto const& cur : all_entries) {
                cur.second.writeFile(out);
            }

            out.close();
        }
    }

    static void addAll(std::map<std::string, DataEntry>& entries, DataEntry const& cur) {
        DataEntry& entry = entries[cur.marker];

        if(0 == entry.calls) {
            entry = cur;
        } else {
            entry.add(cur);
        }
    }
};

extern Perf perf;

struct PerfSentinel {
    std::string marker = "";
    std::chrono::high_resolution_clock::time_point start = {};

    PerfSentinel(char const* file, char const* func, int line) {
        marker += std::string(file) + ":" + std::to_string(line) + ":" + func;
        start = std::chrono::high_resolution_clock::now();
    }

    ~PerfSentinel() {
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;

        perf.add(marker, diff.count());
    }

    PerfSentinel(PerfSentinel const&) = delete;
    PerfSentinel& operator=(PerfSentinel const&) = delete;
};

#define ISSM_PERF_ZONE PerfSentinel _perf_sent(__FILE__, __func__, __LINE__);

#endif  // PERF_HPP
