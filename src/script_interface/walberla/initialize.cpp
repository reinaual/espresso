#include <script_interface/ObjectHandle.hpp>

#include "WalberlaBlockForest.hpp"
#include <utils/Factory.hpp>

namespace ScriptInterface::walberla {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<WalberlaBlockForest>("walberla::WalberlaBlockForest");
}

} // namespace ScriptInterface::walberla