
#ifndef __GOBJECTMANAGER_H__
#define __GOBJECTMANAGER_H__

#include <map>
#include <set>
#include <string>
#include <vector>

#include <TDirectory.h>

class GObjectManager : public TDirectory {
	private:
		GObjectManager();
		static GObjectManager *fGObjectManager;
		std::map<std::string,TObject*> fManagedObjects;
		std::set<std::string> fManagedSources;

	public:
		static GObjectManager *instance();
		~GObjectManager();

		static TObject *Add(TObject* obj);
		static TObject *Get(std::string name);
		static TObject *Get(TObject* obj);
		static TObject *CreateGObject(TObject *obj);
		static TObject *AddObject(const std::string& source,
		                          const std::string& path,
		                          TObject *obj);
		static TObject *UpdateObject(const std::string& source,
		                             const std::string& path,
		                             TObject *obj);
		static TObject *GetObject(const std::string& source,
		                          const std::string& path);
		static bool     HasObject(const std::string& source,
		                          const std::string& path);
		static bool     HasSource(const std::string& source);
		static std::vector<std::string> GetPaths(const std::string& source);
		static std::string MakeKey(const std::string& source,
		                           const std::string& path);
		static unsigned long UpdateCount();


	ClassDef(GObjectManager,0)

};

#endif
