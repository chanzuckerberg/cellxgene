import React from "react";
import { connect } from "react-redux";
import { Button, Tooltip, Dialog, Classes, Colors } from "@blueprintjs/core";

@connect(state => ({
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  ontology: state.ontology,
  ontologyLoading: state.ontology?.loading
}))
class AnnoDialog extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const {
      isActive,
      text,
      title,
      instruction,
      cancelTooltipContent,
      errorMessage,
      validationError,
      annoSelect,
      ontologySelect,
      annoInput,
      handleCancel,
      handleSubmit,
      primaryButtonText,
      secondaryButtonText,
      handleSecondaryButtonSubmit
    } = this.props;

    return (
      <Dialog icon="tag" title={title} isOpen={isActive} onClose={handleCancel}>
        <form
          onSubmit={e => {
            e.preventDefault();
          }}
        >
          <div className={Classes.DIALOG_BODY}>
            <div style={{ marginBottom: 20 }}>
              <p>{instruction}</p>
              {annoInput || null}
              <p
                style={{
                  marginTop: 7,
                  visibility: validationError ? "visible" : "hidden",
                  color: Colors.ORANGE3
                }}
              >
                {errorMessage}
              </p>
            </div>
            {annoSelect || null}
            {ontologySelect || null}
          </div>
          <div className={Classes.DIALOG_FOOTER}>
            <div className={Classes.DIALOG_FOOTER_ACTIONS}>
              <Tooltip content={cancelTooltipContent}>
                <Button onClick={handleCancel}>Cancel</Button>
              </Tooltip>
              {handleSecondaryButtonSubmit && secondaryButtonText ? (
                <Button
                  onClick={handleSecondaryButtonSubmit}
                  disabled={!text || validationError}
                  intent="none"
                  type="submit"
                >
                  {secondaryButtonText}
                </Button>
              ) : null}
              <Button
                onClick={handleSubmit}
                disabled={!text || validationError}
                intent="primary"
                type="submit"
              >
                {primaryButtonText}
              </Button>
            </div>
          </div>
        </form>
      </Dialog>
    );
  }
}

export default AnnoDialog;
