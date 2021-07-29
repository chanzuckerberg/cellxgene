import React from "react";
import { connect } from "react-redux";
import actions from "../../actions";
import FilenameDialog from "./filenameDialog";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  annotations: (state as any).annotations,
  obsAnnotationSaveInProgress:
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (state as any).autosave?.obsAnnotationSaveInProgress ?? false,
  genesetSaveInProgress:
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (state as any).autosave?.genesetSaveInProgress ?? false,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  error: (state as any).autosave?.error,
  writableCategoriesEnabled:
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (state as any).config?.parameters?.annotations ?? false,
  writableGenesetsEnabled: !(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    ((state as any).config?.parameters?.annotations_genesets_readonly ?? true)
  ),
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  annoMatrix: (state as any).annoMatrix,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  genesets: (state as any).genesets,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  lastSavedAnnoMatrix: (state as any).autosave?.lastSavedAnnoMatrix,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  lastSavedGenesets: (state as any).autosave?.lastSavedGenesets,
}))
// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class Autosave extends React.Component<{}, State> {
  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {
      timer: null,
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  componentDidMount() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'writableCategoriesEnabled' does not exis... Remove this comment to see the full error message
    const { writableCategoriesEnabled, writableGenesetsEnabled } = this.props;
    let { timer } = this.state;
    if (timer) clearInterval(timer);
    if (writableCategoriesEnabled || writableGenesetsEnabled) {
      timer = setInterval(this.tick, 2500);
    } else {
      timer = null;
    }
    this.setState({ timer });
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  componentWillUnmount() {
    const { timer } = this.state;
    if (timer) clearInterval(timer);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  tick = () => {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
      dispatch,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsAnnotationSaveInProgress' does not ex... Remove this comment to see the full error message
      obsAnnotationSaveInProgress,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'genesetSaveInProgress' does not exist on... Remove this comment to see the full error message
      genesetSaveInProgress,
    } = this.props;
    if (!obsAnnotationSaveInProgress && this.needToSaveObsAnnotations()) {
      dispatch(actions.saveObsAnnotationsAction());
    }
    if (!genesetSaveInProgress && this.needToSaveGenesets()) {
      dispatch(actions.saveGenesetsAction());
    }
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  needToSaveObsAnnotations = () => {
    /* return true if we need to save obs cell labels, false if we don't */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
    const { annoMatrix, lastSavedAnnoMatrix } = this.props;
    return actions.needToSaveObsAnnotations(annoMatrix, lastSavedAnnoMatrix);
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  needToSaveGenesets = () => {
    /* return true if we need to save gene ses, false if we do not */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'genesets' does not exist on type 'Readon... Remove this comment to see the full error message
    const { genesets, lastSavedGenesets } = this.props;
    return genesets.initialized && genesets.genesets !== lastSavedGenesets;
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  needToSave() {
    return this.needToSaveGenesets() || this.needToSaveObsAnnotations();
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  saveInProgress() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsAnnotationSaveInProgress' does not ex... Remove this comment to see the full error message
    const { obsAnnotationSaveInProgress, genesetSaveInProgress } = this.props;
    return obsAnnotationSaveInProgress || genesetSaveInProgress;
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  statusMessage() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'error' does not exist on type 'Readonly<... Remove this comment to see the full error message
    const { error } = this.props;
    if (error) {
      return `Autosave error: ${error}`;
    }
    return this.needToSave() ? "Unsaved" : "All saved";
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'writableCategoriesEnabled' does not exis... Remove this comment to see the full error message
      writableCategoriesEnabled,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'writableGenesetsEnabled' does not exist ... Remove this comment to see the full error message
      writableGenesetsEnabled,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'lastSavedAnnoMatrix' does not exist on t... Remove this comment to see the full error message
      lastSavedAnnoMatrix,
    } = this.props;
    const initialDataLoadComplete = lastSavedAnnoMatrix;
    if (!writableCategoriesEnabled && !writableGenesetsEnabled) return null;
    return (
      <div
        id="autosave"
        data-testclass={
          !initialDataLoadComplete
            ? "autosave-init"
            : this.saveInProgress() || this.needToSave()
            ? "autosave-incomplete"
            : "autosave-complete"
        }
        style={{
          position: "absolute",
          display: "inherit",
          right: 8,
          bottom: 8,
          zIndex: 1,
        }}
      >
        {this.statusMessage()}
        <FilenameDialog />
      </div>
    );
  }
}

export default Autosave;
