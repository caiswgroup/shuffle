#ifndef _pamap_h_INCLUDED
#define _pamap_h_INCLUDED

#include <iostream>
#include <mutex>
#include <list>
#include <utility>
#include <string>
#include <algorithm>   
#include <memory>
#include "svec.h"
 
 
template<typename Key,typename Value,typename Hash = std::hash<Key>>
class threadsafe_lookup_table
{
private:
	class bucket_type
	{
	public:
		Value value_for(Key const& key,Value const& default_value) const
		{
			_bucket_const_iterator found_entry = find_entry_for(key);
			return (found_entry == _data.end()) ? default_value : found_entry->second;
		}
 
		void add_or_update_mapping(Key const& key,Value const& value)
		{
			_bucket_iterator found_entry = find_entry_for(key);
			if(found_entry == _data.end())
			{
				_data.push_back(_bucket_value(key,value));
			}
			else
			{
				found_entry->second = value;
			}
		}
 
		void remove_mapping(Key const& key)
		{
			_bucket_iterator const found_entry = find_entry_for(key);
			if(found_entry != _data.end())
			{
				_data.erase(found_entry);
			}
		}
	private:
		typedef std::pair<Key,Value> _bucket_value;
		typedef std::list<_bucket_value> _bucket_data;
		typedef typename _bucket_data::iterator _bucket_iterator;
		typedef typename _bucket_data::const_iterator _bucket_const_iterator;
 
		_bucket_data _data;

 
		_bucket_iterator find_entry_for(Key const& key)
		{
			return std::find_if(_data.begin(),_data.end(),
			                [&](_bucket_value const& item)
			                {return item.first==key;});
		}
 
		_bucket_const_iterator find_entry_for(Key const& key) const
		{
			return std::find_if(_data.begin(),_data.end(),
			                [&](_bucket_value const& item)
			                {return item.first==key;});
		}
	};
 
public:
	svec<std::unique_ptr<bucket_type> > buckets;
	Hash hasher;
	bucket_type& get_bucket(Key const& key) const
	{
		std::size_t const bucket_index = hasher(key) % buckets.size();
		return *buckets[bucket_index];
	}
 
public:
	typedef Key key_type;
	typedef Hash hash_type;
	threadsafe_lookup_table(
			unsigned num_buckets = 10000019,Hash const& hasher_ = Hash()):
			buckets(num_buckets),hasher(hasher_)
	{
		for(unsigned i = 0;i < num_buckets;++i)
		{
			buckets[i].reset(new bucket_type);
		}
	}
 
	threadsafe_lookup_table(threadsafe_lookup_table const& other) = delete;
 
	threadsafe_lookup_table const& operator=(
			threadsafe_lookup_table const& other) = delete;
 
	Value value_for(Key const& key,Value const& default_value = Value())
	{
		return get_bucket(key).value_for(key,default_value);
	}
 
	void remove_maping(Key const& key)
	{
		get_bucket(key).remove_mapping(key);
	}
};
 

#endif
